; This reads the NLyc and metallicity from SAG to compute emission lines by interpolating 
; from MAPPING grids
;

pro models,Karr,GammArr,ZstarArr,Nionzmod,modelnameArr
	LineModArr  = [0,0,1,intarr(8)+1]
	BulgeModArr	= [0,0,1,intarr(8)+1]
	qValArr		= ['1e7','4e8','0',strarr(8)+'0']
	betaGMCArr	= ['0','0','0',strarr(8)+'0']
	qVarArr		= ['0','0','2',strarr(8)+'2']
	KArr 		= [1.0e7,4.0e8,9.4e6, 1.6e7, 2.0e7, 7.7e6, 1.3e7, 2.8e7, 5.0e6, 2.2e7, 3.3e7]
	GammArr		= [0,0,-1.0, -0.8, -0.7, -2.0, -1.33, -1.00, -5.3, -2.2, -1.45]
	ZstarArr	= [1,1,0.006, 0.008, 0.014, 0.007, 0.012, 0.013, 0.007, 0.011, 0.017]

	
	KArr 		= [1.0e7,4.0e8, 2.8e7,2.8e7,2.8e7,2.8e7];,2.8e7]
;	GammArr		= [0,0,-0.8, -1.3,-1.8,-2.3]
	GammArr		= [0,0,-1.3,-0.8,-1.8,-2.3]
	ZstarArr	= [1,1,0.012, 0.012,0.012,0.012,0.012]

	modelnameArr = ['q=1e7','q=4e8','\gamma='+strn(-1*GammArr[2],len=3),$
	'\gamma='+strn(-1*GammArr[3],len=3),'\gamma='+strn(-1*GammArr[4],len=3),$
	'\gamma='+strn(-1*GammArr[5],len=3)];,'\gamma='+strn(-1*GammArr[6],len=3)]
	
	Nionzmod = n_elements(KArr)	
;	Nionzmod = 6

	return
end

pro read_mappings,LinesArr,nZ,nQ,nAge,$
	ZArray,QArray,AgeArr,LambdaArr,CentralLambda

	LineDir = '/data3/aaorsi/SAG/Version/sag4_mar/data/Lines/'

;	MappingsFile = LineDir+'LineData_Levesque10_cont_dens_1e1'
;	MappingsInfo = LineDir+'LineInfo_Levesque10_cont_'

	MappingsFile = LineDir+'LineData_Levesque10dens_1e1'
	MappingsInfo = LineDir+'LineInfo_Levesque10'
;	Reading Mappings Data
	print,'Reading Mappings data first'
	
	print,MappingsInfo

	openr,1,MappingsInfo
	str=''
	NumEmLines = numlines(MappingsInfo)-2
	LambdaArr = strarr(NumEmlines)
	CentralLambda = fltarr(NumEmLines)		

	readf,1,str
	readf,1,str
	j = 0
	while ~EOF(1) do begin
		readf,1,str
		res = strsplit(str,/ex)
		LambdaArr[j] = res[1]
		CentralLambda[j] = res[2] + 0.0
		j++
	endwhile
	
	print,'Number of emission lines: ',j,NumEmLines

	print,MappingsFile

	NLines = 0l
	nZ 	   = 0l
	nQ 	   = 0l
	nAge   = 0l

	close,1

	openr,1,MappingsFile
	readu,1,NLines
	readu,1,CentralLambda

	readu,1,nZ
	readu,1,nQ
	readu,1,nAge
	
	ZArray = fltarr(nZ)
	QArray = fltarr(nQ)
	AgeArr = fltarr(nAge)

	readu,1,ZArray
	readu,1,QArray
	readu,1,AgeArr
	
	LinesArr = dblarr(nQ*nZ*nAge*NLines)
	
	readu,1,LinesArr
	
	close,1

	return
end			

; The following is adapted from emlines.c
Function calc_emlines,QArr,ZArr,QGas,ZGas,LinesArray,$
	nZ,nQ,nAgeArr,idl,NLines


	SED_NZEL = nZ
	SED_NQEL = nQ

	age = 0;	
	ia = 0;

	i = 0;

	if (QGas lt 1e7 ) then begin
		line = 1e-30;
;		print,'IonPar = 0, line = 0'
		return,line;
	endif

	if (ZGas gt 1.0 or ZGas lt 1e-4) then begin
		line = 1e-30
;		print,'Z > 1, line = 0'
		return,line
	endif

	if (ZGas le ZArr[0]) then begin
		iz0	= 0;
		iz1 = 1;
		Z0 	= ZArr[iz0];
		Z1	= ZArr[iz1];
	endif else begin
		i = 0;
		while(ZGas ge ZArr[i] and i lt SED_NZEL-1) do begin
			Z0 = ZArr[i];
			iz0 = i;
			i++;
;			if i eq SED_NZEL then break
		endwhile

		iz1 = i 
		Z1 = ZArr[iz1];
	endelse

	if (QGas le QArr[0]) then begin
		iq0	= 0;
		iq1	= 1;
		Q0	= QArr[0];		
		Q1	= QArr[1];
	endif else begin
		i = 0;
		while(QGas ge QArr[i] and i lt SED_NQEL-1)  do begin
			Q0 = QArr[i];	
			iq0 = i;
			i++;
		endwhile
	
		iq1 = i
		Q1 = QArr[iq1];
	endelse

	f00 = LinesArray[idl + NLines*iz0 + SED_NZEL*NLines*iq0];
	f10 = LinesArray[idl + NLines*iz1 + SED_NZEL*NLines*iq0];
	f01 = LinesArray[idl + NLines*iz0 + SED_NZEL*NLines*iq1];
	f11 = LinesArray[idl + NLines*iz1 + SED_NZEL*NLines*iq1];

	den = (Z1-Z0)*(Q1-Q0);
				
	t1 = f00*(Z1-ZGas)*(Q1-QGas)/den;
	t2 = f10*(ZGas-Z0)*(Q1-QGas)/den;
	t3 = f01*(Z1-ZGas)*(QGas-Q0)/den;
	t4 = f11*(ZGas-Z0)*(QGas-Q0)/den;

	line = t1 + t2 + t3 + t4;
	line = (line le 0.0) ? 1e-30 : line;

	return,line
end

pro integ_line,QArr,ZArr,QGas,ZGas,LinesArr,$
	nZ,nQ,nAgeArr,nlyc,idl,NLines,line,lineSAG=lineSAG

	QMin = 1.0e7
	QMax = 4.0e8
	
; Q-Z Model

;	K		= Kqz;
;	M_		= Zqz;
;	Gamm	= Gqz;
		
;	Mean_q = K * pow(ZGal/M_,Gamm);

	Mean_q = QGas

;	if (Mean_q lt QMin) then Mean_q = QMin;
;	if (Mean_q gt QMax) then Mean_q = QMax;
	
;	ZGas *= pow(10,Zoff);

	FracHyda   = calc_emlines(QArr,ZArr,Mean_q,ZGas,LinesArr,$
	nZ,nQ,nAgeArr,2,NLines)

	cel = calc_emlines(QArr,ZArr,Mean_q,ZGas,LinesArr,$
	nZ,nQ,nAgeArr,idl,NLines)

	alpha = alog10(1.37) - 12.0
	line = (cel le 1e-29 or FracHyda le 1e-29 or ~Finite(nlyc)) ? -30.00 :  $
			nlyc + alpha + alog10(cel) - alog10(FracHyda) -40.0;
	
	if keyword_set(LineSag) then $
		if (abs(10^line - lineSAG)/lineSAG gt 1e-3) then stop

	if nlyc lt 20 and line gt -2 then stop,'Nlyc too low, line too high'

	return


end
	
pro get_emlines
	print,'Emission lines procedures compiled'
	return
end

