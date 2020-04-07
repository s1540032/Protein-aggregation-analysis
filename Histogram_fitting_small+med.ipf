#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.



function pathtosavedata()
setdatafolder root:
make/o/n=1/T pathway="Macintosh HD:Users:Mathew:Desktop:Output:"		/// Change output folder

end

macro runmedium()
pathtosavedata()
loadmedium()
trialfitmedium()
final_fit_medium()
saveallmedium()

end

function loadmedium() 
setdatafolder root:
newdatafolder/s medium
// Function to load the FRET data into one table
LoadWave/J/M/D/A=wave/K=0 "Clipboard"

// This is loaded into a 2D matrix


end

function trialfitmedium()
wave wave0		// Let knows exist. 

make/o/n=(dimsize(wave0,0)) FRET_axis			// Need to make the FRET axis for the plot

variable a,b,c,d									// Some variables

for(a=0;a<(dimsize(FRET_axis,0));a+=1)		// Go through each row
	FRET_axis[a]=wave0[a][0]					// Populate from the first column
endfor



variable number_of_fits=(dimsize(wave0,1)-1) 		// Number of time-points to fit. 
variable length=(dimsize(wave0,0))
make/o/n=(length) tempy			// Temp wave to store data in


make/o/n=(number_of_fits) y0,a1,x1=NAN,w1=NAN,a2,x2=NAN,w2=NAN		// To store variables in
make/o/n=1 x1_ave,x2_ave,w1_ave,w2_ave		//

for(a=0;a<(number_of_fits);a+=1)				// Go through and fit them all

	//First of all need to extract the data
	for(b=0;b<length;b+=1)
	
	tempy[b]=wave0[b][a+1]
	
	
	endfor
	
	string name="Y_"+num2str(a)			// Duplicate the wave to store
	duplicate/o tempy,$name
	make/o/T T_Constraints= {"K1 > 0","K2 > 0.05","K2 < 0.5","K3 > 0.05","K3 < 0.5","K4 > 0","K5 > 0.5","K6 > 0.05","K6 < 0.5"}
	make/o/n=7 W_coef= {0,10,0.5,0.1,10,0.8,0.1}		// Perform the fits on all of the data
	FuncFit twogauss W_coef tempy /X=FRET_axis /D/C=T_Constraints
	
	// Store the fits in a table.
	
	y0[a]=w_coef[0]
	a1[a]=w_coef[1]
	if(w_coef[2]<0.5 && w_coef[2]>0.05)
	x1[a]=w_coef[2]
	endif
	if(w_coef[3]<0.5 && w_coef[3]>0.05)
	w1[a]=w_coef[3]
	endif
	a2[a]=w_coef[4]
	if(w_coef[5]>0.5 && w_coef[5]<1)
	x2[a]=w_coef[5]
	endif
	if(w_coef[6]<0.5 && w_coef[6]>0.05)
	w2[a]=w_coef[6]
	endif

endfor

// Now take an average for the next fits

wavestats/q x1
x1_ave[0]=v_avg
wavestats/q w1
w1_ave[0]=v_avg



wavestats/q x2
x2_ave[0]=v_avg
wavestats/q w2
w2_ave[0]=v_avg


end



function final_fit_medium()
wave wave0		// Let knows exist. 

variable a,b,c,d									// Some variables

make/o/n=200 fit_xaxis // This is for the fit axis
variable e=0
for(a=0;a<200;a+=1)
fit_xaxis[a]=e
e+=0.005
endfor



wave FRET_axis




variable number_of_fits=(dimsize(wave0,1)-1) 		// Number of time-points to fit. 
variable length=(dimsize(wave0,0))
make/o/n=(length) tempy			// Temp wave to store data in


make/o/n=(number_of_fits) y0_fit,a1_fit,x1_fit,w1_fit,a2_fit,x2_fit,w2_fit,int1,int2,conc1,conc2		// To store variables in

wave x1_ave,x2_ave,w2_ave,w1_ave		// This is to get the data for the fits. 

for(a=0;a<(number_of_fits);a+=1)				// Go through and fit them all

	//First of all need to extract the data
	for(b=0;b<length;b+=1)
	
	tempy[b]=wave0[b][a+1]
	
	
	endfor
	
	string name="Y_"+num2str(a)			// Duplicate the wave to store
	duplicate/o tempy,$name

	make/o/n=7 W_coef= {0,10,0.3,0.1,10,0.5,0.1}		// Perform the fits on all of the data
	w_coef[0]=0
	w_coef[1]=10
	w_coef[2]=x1_ave[0]
	w_coef[3]=w1_ave[0]
	w_coef[4]=10
	w_coef[5]=x2_ave[0]
	w_coef[6]=w2_ave[0]
	
	Make/O/T/N=2 T_Constraints = {"K1 > 0","K4 > 0"}
	FuncFit/H="1011011" twogauss W_coef tempy /X=FRET_axis /D/C=T_Constraints 
	
	string fitsave="total_fit_"+num2str(a)
	wave fit_tempy
	duplicate/o fit_tempy,$fitsave
	// Store the fits in a table.
	
	y0_fit[a]=w_coef[0]
	a1_fit[a]=w_coef[1]
	x1_fit[a]=w_coef[2]
	w1_fit[a]=w_coef[3]
	a2_fit[a]=w_coef[4]
	x2_fit[a]=w_coef[5]
	w2_fit[a]=w_coef[6]
	
	// Now for single gauss fits to add to graphs:
	make/o/n=4 w_coef
	
	k0=y0_fit[a]
	k1=a1_fit[a]
	k2=x1_fit[a]
	k3=w1_fit[a]
	
	CurveFit/H="1111" gauss tempy /X=FRET_axis /D 
	
	string peak1="Peak_1_"+num2str(a)
	duplicate/o fit_tempy,$peak1
	
	k0=y0_fit[a]
	k1=a2_fit[a]
	k2=x2_fit[a]
	k3=w2_fit[a]
	
	CurveFit/H="1111" gauss tempy /X=FRET_axis /D 
	
	string peak2="Peak_2_"+num2str(a)
	duplicate/o fit_tempy,$peak2
	
	Integrate $peak1/D=peak_INT;DelayUpdate
	
	int1[a]=peak_int[199]
	conc1[a]=int1[a]/0.05
	
	Integrate $peak2/D=peak2_INT;DelayUpdate
	
	int2[a]=peak2_int[199]
	conc2[a]=int2[a]/0.05
	
	
	
	display $name vs FRET_axis	
	ModifyGraph mode=5,rgb=(39321,39321,39321)
	ModifyGraph width=283.465,height=198.425
	ModifyGraph gFont="Arial",gfSize=20
	Label left "No. of oligomers";DelayUpdate
	Label bottom "FRET Efficiency"
	ModifyGraph height=170.079
	AppendToGraph $fitsave vs fit_xaxis
	ModifyGraph lsize($fitsave)=2,rgb($fitsave)=(1,3,39321)
	AppendToGraph $peak1,$peak2 vs fit_xaxis
	ModifyGraph lsize($peak2)=2,rgb($peak2)=(65535,43690,0),lsize($peak1)=2;DelayUpdate
	ModifyGraph rgb($peak1)=(2,39321,1)
	
	
endfor




end


function see(number,maxsmall,maxmed)
variable number,maxsmall,maxmed
setdatafolder root:
setdatafolder medium
wave fret_axis,fit_xaxis
string name="Y_"+num2str(number)
string fitsave="total_fit_"+num2str(number)
string peak1="Peak_1_"+num2str(number)
string peak2="Peak_2_"+num2str(number)



display $name vs FRET_axis	
	ModifyGraph mode=5,rgb=(39321,39321,39321)
	ModifyGraph width=283.465,height=198.425
	ModifyGraph gFont="Arial",gfSize=20
	Label left "No. of oligomers";DelayUpdate
	Label bottom "FRET Efficiency"
	ModifyGraph height=170.079
	AppendToGraph $fitsave vs fit_xaxis
	ModifyGraph lsize($fitsave)=2,rgb($fitsave)=(1,3,39321)
	AppendToGraph $peak1,$peak2 vs fit_xaxis
	ModifyGraph lsize($peak2)=2,rgb($peak2)=(65535,43690,0),lsize($peak1)=2;DelayUpdate
	ModifyGraph rgb($peak1)=(2,39321,1)
	SetAxis left *,maxmed
setdatafolder root:
setdatafolder small

display $name vs FRET_axis	
	ModifyGraph mode=5,rgb=(39321,39321,39321)
	ModifyGraph width=283.465,height=198.425
	ModifyGraph gFont="Arial",gfSize=20
	Label left "No. of oligomers";DelayUpdate
	Label bottom "FRET Efficiency"
	ModifyGraph height=170.079
	SetAxis left *,maxsmall

	AppendToGraph $peak1 vs fit_xaxis
	ModifyGraph lsize($peak1)=2;DelayUpdate
	ModifyGraph rgb($peak1)=(1,3,39321)


end



function saveallmedium()
setdatafolder root:
wave/t pathway
string path=pathway[0]
setdatafolder medium

wave wave0
variable a

for(a=0;a<(dimsize(wave0,1)-1);a+=1)
variable number=a
wave fret_axis,fit_xaxis
string name="Y_"+num2str(number)
string fitsave="total_fit_"+num2str(number)
string peak1="Peak_1_"+num2str(number)
string peak2="Peak_2_"+num2str(number)



display $name vs FRET_axis	
	ModifyGraph mode=5,rgb=(39321,39321,39321)
	ModifyGraph width=283.465,height=198.425
	ModifyGraph gFont="Arial",gfSize=20
	Label left "No. of oligomers";DelayUpdate
	Label bottom "FRET Efficiency"
	ModifyGraph height=170.079
	AppendToGraph $fitsave vs fit_xaxis
	ModifyGraph lsize($fitsave)=2,rgb($fitsave)=(1,3,39321)
	AppendToGraph $peak1,$peak2 vs fit_xaxis
	ModifyGraph lsize($peak2)=2,rgb($peak2)=(65535,43690,0),lsize($peak1)=2;DelayUpdate
	ModifyGraph rgb($peak1)=(2,39321,1)
	
	string path2=path+"Histogram_Small_"+num2str(a)+".pdf"
	SavePICT/E=-8/EF=1/o as path2
	endfor

Edit/K=0 root:medium:conc1
Edit/K=0 root:medium:conc2
end






Function twogauss(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0 + a1*exp(-((x-xc1)/w1)^2) + a2*exp(-((x-xc2)/w2)^2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 7
	//CurveFitDialog/ w[0] = y0
	//CurveFitDialog/ w[1] = a1
	//CurveFitDialog/ w[2] = xc1
	//CurveFitDialog/ w[3] = w1
	//CurveFitDialog/ w[4] = a2
	//CurveFitDialog/ w[5] = xc2
	//CurveFitDialog/ w[6] = w2

	return w[0] + w[1]*exp(-((x-w[2])/w[3])^2) + w[4]*exp(-((x-w[5])/w[6])^2)
End


function kill()

	variable                      winMask;
 
	variable                      i,n;
	variable                      all=0x1000+0x40+0x10+0x4+0x2+0x1;
	string                        theWins;
 
	winMask = !winMask ? all : winMask;
 
	theWins = winList("*",";","WIN:"+num2iStr(winMask & all));
	for(i=0,n=itemsInList(theWins,";") ; i<n ; i+=1)
		doWindow/K $stringFromList(i,theWins,";");
	endfor;
end








///////////////////SMALL//////////////////////////

macro runsmall()
loadsmall()
trialfitsmall()
final_fit_small()
saveallsmall()

end

function loadsmall() 
setdatafolder root:
newdatafolder/s small
// Function to load the FRET data into one table
LoadWave/J/M/D/A=wave/K=0 "Clipboard"

// This is loaded into a 2D matrix


end

function trialfitsmall()
wave wave0		// Let knows exist. 

make/o/n=(dimsize(wave0,0)) FRET_axis			// Need to make the FRET axis for the plot

variable a,b,c,d									// Some variables

for(a=0;a<(dimsize(FRET_axis,0));a+=1)		// Go through each row
	FRET_axis[a]=wave0[a][0]					// Populate from the first column
endfor



variable number_of_fits=(dimsize(wave0,1)-1) 		// Number of time-points to fit. 
variable length=(dimsize(wave0,0))
make/o/n=(length) tempy			// Temp wave to store data in


make/o/n=(number_of_fits) y0,a1,x1=NAN,w1=NAN
make/o/n=1 x1_ave,w1_ave		//

for(a=0;a<(number_of_fits);a+=1)				// Go through and fit them all

	//First of all need to extract the data
	for(b=0;b<length;b+=1)
	
	tempy[b]=wave0[b][a+1]
	
	
	endfor
	
	string name="Y_"+num2str(a)			// Duplicate the wave to store
	
	
	make/o/n=4 W_coef= {0,10,0.5,0.1}		// Perform the fits on all of the data
	
	
	CurveFit/H="1000" gauss kwCWave=W_coef, tempy /X=FRET_axis /D 
	
	// Store the fits in a table.
	
	y0[a]=w_coef[0]
	a1[a]=w_coef[1]
	x1[a]=w_coef[2]
	w1[a]=w_coef[3]


endfor

// Now take an average for the next fits

wavestats/q x1
x1_ave[0]=v_avg
wavestats/q w1
w1_ave[0]=v_avg


end



function final_fit_small()
wave wave0		// Let knows exist. 

variable a,b,c,d									// Some variables

make/o/n=200 fit_xaxis // This is for the fit axis
variable e=0
for(a=0;a<200;a+=1)
fit_xaxis[a]=e
e+=0.005
endfor



wave FRET_axis




variable number_of_fits=(dimsize(wave0,1)-1) 		// Number of time-points to fit. 
variable length=(dimsize(wave0,0))
make/o/n=(length) tempy			// Temp wave to store data in


make/o/n=(number_of_fits) y0_fit,a1_fit,x1_fit,w1_fit,int1,conc1	// To store variables in

wave x1_ave,w1_ave		// This is to get the data for the fits. 

for(a=0;a<(number_of_fits);a+=1)				// Go through and fit them all

	//First of all need to extract the data
	for(b=0;b<length;b+=1)
	
	tempy[b]=wave0[b][a+1]
	
	
	endfor
	
	string name="Y_"+num2str(a)			// Duplicate the wave to store
	duplicate/o tempy,$name

	make/o/n=4 W_coef= {0,10,0.5,0.1}		// Perform the fits on all of the data
	w_coef[0]=0
	w_coef[2]=x1_ave[0]
	w_coef[3]=w1_ave[0]
	
	CurveFit/H="1011" gauss kwCWave=W_coef, tempy /X=FRET_axis /D 
	
	
	
	string fitsave="total_fit_"+num2str(a)
	wave fit_tempy
	duplicate/o fit_tempy,$fitsave
	// Store the fits in a table.
	
	y0_fit[a]=w_coef[0]
	a1_fit[a]=w_coef[1]
	x1_fit[a]=w_coef[2]
	w1_fit[a]=w_coef[3]
	
	

	
	string peak1="Peak_1_"+num2str(a)
	duplicate/o fit_tempy,$peak1
	
	
	Integrate $peak1/D=peak_INT;DelayUpdate
	
	int1[a]=peak_int[199]
	conc1[a]=int1[a]/0.05
	

	
	
	display $name vs FRET_axis	
	ModifyGraph mode=5,rgb=(39321,39321,39321)
	ModifyGraph width=283.465,height=198.425
	ModifyGraph gFont="Arial",gfSize=20
	Label left "No. of oligomers";DelayUpdate
	Label bottom "FRET Efficiency"
	ModifyGraph height=170.079

	AppendToGraph $peak1 vs fit_xaxis
	ModifyGraph lsize($peak1)=2;DelayUpdate
	ModifyGraph rgb($peak1)=(1,3,39321)
	
	
endfor




end




function saveallsmall()
setdatafolder root:
wave/t pathway
string path=pathway[0]
setdatafolder small

wave wave0
variable a

for(a=0;a<(dimsize(wave0,1)-1);a+=1)
variable number=a
wave fret_axis,fit_xaxis
string name="Y_"+num2str(number)

string peak1="Peak_1_"+num2str(number)




display $name vs FRET_axis	
	ModifyGraph mode=5,rgb=(39321,39321,39321)
	ModifyGraph width=283.465,height=198.425
	ModifyGraph gFont="Arial",gfSize=20
	Label left "No. of oligomers";DelayUpdate
	Label bottom "FRET Efficiency"
	ModifyGraph height=170.079

	AppendToGraph $peak1 vs fit_xaxis
	ModifyGraph lsize($peak1)=2;DelayUpdate
	ModifyGraph rgb($peak1)=(1,3,39321)
	
	string path2=path+"Histogram_Small"+num2str(a)+".pdf"
	SavePICT/E=-8/EF=1/o as path2
	endfor

Edit/K=0 root:small:conc1

end



macro close_windows()
kill()
end



Macro Look(Timepoint,small,med)
variable small=1000
Prompt Small, "Smallmers axis max: "
variable med=1000
Prompt med, "Mediummers axis max: "
variable Timepoint=0
Prompt Timepoint, "Timepoint to see: "			// Set prompt for X param


see(timepoint,small,med)

end