//To be used in conjunction with FPGA card and output .DAT files. Copyright Mathew Horrocks, 2012. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//List of Commands:
//	loader()										Change variables- loads files from FPGA card.
//	experiment(small,medium,large)				Generates tables used to analyse the data. 
//	kill()											Kills windows for display function. 
//	ana(time)									Master command for each time-point to be analysed.
//	get()										Loads data and calculates FRET efficiency and sizes.
// 	matrices_fret()								Generates 2D FRET histograms
//	fret_histograms()								Generates FRET histograms for each different size. 
// 	matrices_lnz()								Generates 2D Lnz histograms
// 	see(path)									Generates print-out. 
//	singlehist(path)								Generates FRET histogram for each size. 
//	print2d(z)									Layout all 2D number distributions on one page. z is number of time-points. 
//	totalsfret()									Generates graph of total populations over time. 
// 	twod()										Layout all 2D number distributions. 
//	mass()										Layout all mass distributions. 
//	error()										Deletes last time-point- must also manually delete the folder. 
//	count()										Counter to keep track. 
//	donorbursts()									Monomer counts
//	Sizes()										Graph showing size-distributions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma rtGlobals=1		// Use modern global access method.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Experiment initiator/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// This function initiates the program.

Function experiment(file,donor,acceptor,first,last,bins,sma,med,binwidth)
variable sma,med,donor,acceptor,first,last,bins,binwidth		// Variables loaded via macro. 

string file

///////// Files to load here //////////

variable numberofexperiments=1							// Change here for the number of experiments
variable exps=numberofexperiments
make/o/T/n=(numberofexperiments) filelist	


filelist[0]="SARAH:Confocal:Sarah:20200225:6hr_1_20nM_FRET:"


variable m,n
//killwaves/a												// Kill all waves
variable b = 20      										// BIN SIZE OVER WHICH DATA SHOULD BE REMOVED

make/o/n=(numberofexperiments) timepoints									// Timepoints made. 	
for(m=0;m<(numberofexperiments);m+=1)
	timepoints[m]=m
endfor
make/o/n=(dimsize(wave0,0)) Width_deletion					// Makes a folder to store deleted events. 
duplicate wave1,times2									// Times. 
killwaves wave0											// Kill waves- keep tide. 								
killwaves wave1
wave timepoints
variable x
variable names
string name
						// Tells program how long to make the tables. 
make/o/n=5 Variables
Variables[1] = sma										// Small cutoff boundary
Variables[2] = med										// Medium cutoff boundary. 

// Make tables to store data. 
	
Make/o/n=(60,exps+1) lnz_Small, lnz_Medium, lnz_Large
Make/o/n=(20,exps+1) FRET_Small, FRET_Medium, FRET_Large,AllFret,verylarge
make/o/n=1 counter=0

// Make FRET axis scale. 

variable r = 0.025
	for(m=0;m<=19;m+=1)
		FRET_Small[m][0] = r
		FRET_Medium[m][0] = r
		FRET_Large[m][0] = r
		VeryLarge[m][0]=r
		allfret[m][0]=r
		r+=0.05
	endfor

// Make z-parameter axis-scale. 

variable t = -2.95
	for(m=0;m<=60;m+=1)
		lnz_Small[m][0] = t
		lnz_Medium[m][0] = t
		lnz_Large[m][0] = t
		t+=0.1
	endfor
variable pod
wave timepoints
// This cycles through the program for each listing in the index file. 	
	do
		name=num2str(x)		// Gets hours from the index. 
		names =x			// Gets hours as variable
		string pathname=filelist[x]
		count(names)					// Adds count to counter
		string nameit=num2str(names)		

		newdatafolder/s/o  $name		// Make new data folder with name as hours- set as current folder

		copy()						// Copy necessary waves into this folder. 	

		make/o/n=1 hour=names		// Need hours within folder
		variables[4]=names
		wave variables				
		// Makes variables which are needed for future functions. 
		make/o/n=5 Variables
		Variables[1] = sma
		Variables[2] = med
		variables[4]=names
		string pathway=pathname
		//string pathway=path+pathname
		loader(pathname,file,donor,acceptor,first,last,bins,binwidth) // Loader
		get()		// Extracts data. 
		consecu()	// Deletes events which are too large. 
		zplot()		// Makes z-parameter histograms. 
		matrices_fret()	// Makes 2D plots
		FRET_histograms()	// Makes FRET efficiency histograms. 
		matrices_z()
		lnz_histograms()
		totalsFRET()			// Counts events. 
		x+=1
		pod+=1
	while(x<numberofexperiments)

totalsfret()

end




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Loader for files ////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function loads the files into the program- works with .DAT files from FPGA cards on instruments. 
// Change variables here- i.e. gamma factor etc.

function loader(pathway,filename,donorthresh,acceptorthresh,first,last,bins,binwidth)
// These are variables and strings used by the program- they are loaded from experiment()
string pathway 				// Folder pathway. 
variable donorthresh			// Donor threshold
variable acceptorthresh		// Acceptor threshold
variable first					// First filenumber
variable last					// Last filenumber
variable bins					// Number of bins per file
variable binwidth				// Bin-width
string filename				// File-name
// CHANGE THESE FOR EACH INSTRUMENT:
variable autodonor=0.200		// Autofluorescence in the donor channel- i.e. average intensity from buffer.
variable autoacceptor=0.15	// Autofluorescence in acceptor channel. 
variable crosstalk=0.107		// Crosstalk- i.e. some donor intensity passes into acceptor channel- this is decimal
variable gammafactor=1.0		// Gamma factor to correct for detection efficiencies and quantum yields. 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
variable m,n=0,g,h
variable d=bins*(last-first) 		// Number of bins in total
variable frametime=(binwidth/1000000000)*bins	// Total frame time
make/o/n=1 gamma_factor=gammafactor		// Make folder with gamma factor in it so that other functions can refer back to it. 
make/o/n=(1) FRET_donor=NAN,FRET_acceptor=NAN,donorchannel=NAN,bin_number=NAN		// Makes waves to store events in. 
make/o/n=((last-first),3) Stats	// Wave for statistics. 

// This section is to get the correct filename and pathway to the files to be loaded. 


string num
variable c



for(c=0;c<last;c+=1)
	variable ee=c+1
	if(ee<2)
		string nametoload=filename
	elseif(ee<10)
		nametoload=filename+"_0"+num2str(ee)
	else
		nametoload=filename+"_"+num2str(ee)
	endif

	// Following to load files

	string namer=num2str(c)
	string loader=pathway+nametoload
	printf "File to load is %s \r",loader
	LoadWave/J/D/W/K=0/A loader





	wave wave0,wave1,wave2				// Tell program that files have been loaded, can then refer back to them.
	variable e,k,f
		for(e=0;e<=(dimsize(wave0,0));e+=1)	// Go through all of bins
			if(wave0[e]>donorthresh&&acceptorthresh<(wave1[e]-crosstalk*wave0[e])) 	// Thresholding- only select files greater than threshold.
				Redimension/N=(k) FRET_Donor,FRET_Acceptor,bin_number			// Increase size of waves by 1 to add data.
				FRET_Donor[k]=wave0[e]-autodonor								// Add data. 
				FRET_Acceptor[k]=wave1[e]-autoacceptor-crosstalk*wave0[e]
				bin_number[k]=e
				k+=1
				m+=1
			elseif(wave0[e]>donorthresh)					// This part takes those bursts that are only in the donor channel- monomeric bursts- needed for the "monomer brightness"
				Redimension/N=(f) donorchannel
				donorchannel[f]=wave0[e]
				f+=1
				n+=1
			endif
		endfor
	stats[c][0]=n/frametime	// Statistics.
	stats[c][1]=m/frametime
	stats[c][2]=m/n
	m=0
	n=0
	killwaves wave0			// Kill waves to allow for load of the next frame and selection of events. 
	killwaves wave1
	killwaves wave2

endfor

display stats[][0]				// Shows the statistics (donor bursts)
// Following makes graph look nice:
TextBox/C/N=text0/F=0/A=MT/E Pathway	
Label bottom "Frame Number"
Label left "Event rate (s\\S-1\\M)"
AppendToGraph/R stats[][1]
Label right "Event rate (s\\S-1\\M)"
ModifyGraph rgb(Stats#1)=(0,0,65280)
ModifyGraph axRGB(right)=(0,0,65280)
ModifyGraph axRGB(left)=(65280,3328,0),tlblRGB(left)=(65280,3328,0);DelayUpdate
ModifyGraph tlblRGB(right)=(0,0,65280),alblRGB(left)=(65280,3328,0);DelayUpdate
ModifyGraph alblRGB(right)=(0,0,65280)
make/o/n=1 donor_bursts=k+f
wavestats donorchannel
make/o/n=1 monomerb=V_avg
SetAxis left 0,*
SetAxis right 0,*
Legend/C/N=text1/J/F=2/A=RC/E "\\s(Stats) Donor\r\\s(Stats#1) Coincident"
print V_avg
//
end









////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////// Kill windows function ///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to close all windows- needed for preview. 
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







//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////Copier///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 // Copies necessary waves from root into new folders. 
 
function copy()
wave width_deletion,variables, FRET_Small,FRET_Small_mass,FRET_medium_mass,FRET_large_mass,FRET_medium,FRET_large,Lnz_small,lnz_medium,lnz_large,counter
duplicate/o root:variables,Variables
duplicate/o root:FRET_Small,FRET_Small
duplicate/o root:FRET_Medium, FRET_Medium
duplicate/o root:FRET_Large, FRET_Large
duplicate/o root:lnz_Small, lnz_Small
duplicate/o root:lnz_Medium, lnz_Medium
duplicate/o root:lnz_Large, lnz_Large
duplicate/o root:counter, counter1
duplicate/o root:width_deletion, width_deletion1
duplicate/o root:verylarge, verylarge
duplicate/o root:allfret, allfret
end

// Copies files when using add function. 
function copy2()
wave width_deletion,variables, FRET_Small,FRET_Small_mass,FRET_medium_mass,FRET_large_mass,FRET_medium,FRET_large,Lnz_small,lnz_medium,lnz_large,counter
duplicate/o root:variables,Variables
duplicate/o root:FRET_small, FRET_small
duplicate/o root:FRET_Medium, FRET_Medium
duplicate/o root:FRET_Large, FRET_Large
duplicate/o root:lnz_small, lnz_small
duplicate/o root:lnz_Medium, lnz_Medium
duplicate/o root:lnz_Large, lnz_Large
end




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////Extracts data from files/////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function get()
variable j
wave variables
wave monomerb
variable mb=monomerb[0] // Monomer brightness. 
wave gamma_factor
variable gammafact=gamma_factor[0]
make/o/n=1 mb1=mb
wave FRET_Acceptor,FRET_Donor,FRET_Efficiency,File_number,Bin_number // Loads waves needed for function- generated from loader()
//Make tables for use later
Make/O/N=(Dimsize(FRET_donor,0),9) Data
Make/O/N=(Dimsize(FRET_donor,0),1) FRET_Efficiencies

//Copy data to another table, and calculate the size and FRET Efficiency
	for(j=0;j<=(Dimsize(FRET_donor,0));j+=1)
	                      
		Data[j][0] = j	  // Width of Bin	
		Data[j][1] = FRET_Donor[j]  // Donor channel intensity
		Data[j][2] = FRET_Acceptor[j]  // Acceptor channel intensity
		Data[j][8] = ln(FRET_Acceptor[j]/FRET_Donor[j]) //ln(Ia/Id)
		Data[j][5] = bin_number[j]  // bin number
		Data[j][6] = fret_acceptor[j]/(gammafact*fret_donor[j]+fret_acceptor[j]) // FRET Efficiencies
		Data[j][7] = (2*((FRET_Donor[j]+(FRET_Acceptor[j]/gammafact))/mb))  // Sizes
	
	Endfor 
Make/N=30/O FRET_Efficiencies_Hist;DelayUpdate
make/N=20/O FRET_Efficiencies_Hist;DelayUpdate
Histogram/C/B={0,0.05,20} FRET_Efficiencies,FRET_Efficiencies_Hist
end





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////Make 2D FRET Histograms///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function matrices_fret()
wave data,A_bursts
variable l,count1,count2,k,count5,m
// Code to make 2D matrices
//Make the appropriate tables for adding data to
Make/O/N=(20,30) Matrix = 0
Make/O/N=(Dimsize(data,0),8) FretH

//Bin the size and FRET data
for(l=0;l<=(Dimsize(data,0));l+=1)
	FretH[l][0] = data[l][6]/0.05				// FRET Efficiency
	FretH[l][1] = data[l][7]/5                 		// Size
	FretH[l][2] = round(FretH[l][0])			// Rounded FRET Efficiency
	FretH[l][3] = round(FretH[l][1])			// Rounded Size
	FretH[l][4] = FretH[l][0] - FretH[l][2]		// Needed to convert rounded to int
	FretH[l][5] = FretH[l][1] - FretH[l][3]		// Needed to convert rounded to int
if(FretH[l][4] < 0)
	FretH[l][6] = FretH[l][2] -1				// FRET Efficiency Int
else
	FretH[l][6] = FretH[l][2]				// FRET Efficiency Int
endif
if(FretH[l][5] < 0)
	FretH[l][7] = FretH[l][3] -1				// Size Int
else
	FretH[l][7] = FretH[l][3]				// Size Int
endif
endfor

////The following is data not accounting for Pascal

//Delete any data in which the size doesn't fit into the matrix
for(k=0;k<=(Dimsize(FretH,0));k+=1)
		if(FretH[k][6] > 20)
		DeletePoints k,1, FretH
		count1 +=1
		endif
endfor
for(k=0;k<=(Dimsize(FretH,0));k+=1)
	if(FretH[k][7] > 30)
	DeletePoints k,1, FretH
	count2 +=1
	endif
endfor
Printf "%d data points deleted due to FRET being too high.\r",count1
Printf "%d data points deleted due to Size being too high.\r",count2
make/o/n=1 toobig=count2

//Add the data to the matrix.
variable x,y
for(m=0;m<=(Dimsize(FretH,0));m+=1)
if(FretH[m][7]<=29)
if(FretH[m][6]<=19)
y = FretH[m][6]
x = FretH[m][7] 
duplicate/o matrix, matrix10
matrix[y][x]=matrix10[y][x]+1 // Adds 1 each time data needs to be added
endif
endif
endfor

////Accounting for Pascal
Make/O/N=(20,150) pascalm = 0
Make/O/N=(Dimsize(data,0),8) Pascal_Fret
for(l=0;l<=(Dimsize(data,0));l+=1)
	Pascal_Fret[l][0] = data[l][6]/0.05				// Bin FRET Efficiency
	Pascal_Fret[l][1] = data[l][7]					//sizes not binned. 
	Pascal_Fret[l][2] = round(Pascal_Fret[l][0])
	Pascal_Fret[l][3] = round(Pascal_Fret[l][1])
	Pascal_Fret[l][4] = Pascal_Fret[l][0] - Pascal_Fret[l][2]
	Pascal_Fret[l][5] = Pascal_Fret[l][1] - Pascal_Fret[l][3]
if(Pascal_Fret[l][4] < 0)
	Pascal_Fret[l][6] = Pascal_Fret[l][2] -1
else
	Pascal_Fret[l][6] = Pascal_Fret[l][2]		// FRET Efficiency
endif
if(Pascal_FRET[l][5] < 0)
	Pascal_Fret[l][7] = Pascal_Fret[l][3] -1
else
	Pascal_Fret[l][7] = Pascal_Fret[l][3]		// Sizes
endif
endfor

//Delete points that won't fit into the 150x20 Matrix

for(k=0;k<=(Dimsize(Pascal_Fret,0));k+=1)
		if(Pascal_Fret[k][6] > 20)
		DeletePoints k,1, Pascal_Fret
		endif
endfor

for(k=0;k<=(Dimsize(Pascal_Fret,0));k+=1)
	if(pascal_FRET[k][7] > 150)
	DeletePoints k,1, Pascal_Fret
	endif

endfor

//Put data into the 150x20 Matrix

for(m=0;m<=(Dimsize(Pascal_Fret,0));m+=1)
variable count4
if(pascal_FRET[m][7]<=149)
if(pascal_FRET[m][6]<=19)
y = Pascal_Fret[m][6]
x = Pascal_Fret[m][7] 
duplicate/o pascalm, pascalm10
pascalm[y][x]=pascalm10[y][x]+1
endif
endif
endfor


//Multiply size distributions by Pascal factor to account for "invisible" oligomers

for(m=0;m<=(Dimsize(pascalm,0));m+=1)

// This accounts for pascal- I have commented out for the time being, since our size estimation is quite poor!

//pascalm[m][2]=pascalm[m][2]*2
//pascalm[m][3]=pascalm[m][3]*1.33
//pascalm[m][4]=pascalm[m][4]*1.14
//pascalm[m][5]=pascalm[m][5]*1.067



endfor


////Now to modify the matrix not accounted for Pascal
duplicate/o matrix, Matrix_pascal
for(m=0;m<=(Dimsize(Matrix_pascal,0));m+=1)

matrix_pascal[m][0]=pascalm[m][0]+pascalm[m][1]+pascalm[m][2]+pascalm[m][3]+pascalm[m][4]
matrix_pascal[m][1]=pascalm[m][5]+pascalm[m][6]+pascalm[m][7]+pascalm[m][8]+pascalm[m][9]

endfor
duplicate/o pascalm FRET_Matrix_1Bin
killwaves pascalm
killwaves pascalm10
//killwaves Pascal_FRET
killwaves matrix10
killwaves FretH
killwaves matrix
duplicate matrix_pascal FRET_Matrix




duplicate/o FRET_Matrix_1Bin FRET_Matrix_1Bin_Mass
variable n=1
for(k=0;k<=(Dimsize(FRET_Matrix_1Bin,1));k+=1)
for(m=0;m<=(Dimsize(FRET_Matrix_1Bin,0));m+=1)

FRET_Matrix_1Bin_Mass[m][k] = FRET_Matrix_1Bin[m][k]*k

n+=1
endfor
endfor
make/o/n=20 FRETaxis
make/o/n=150 size1axis
variable r = 0.025
variable t = 0
for(m=0;m<=19;m+=1)

FRETaxis[m] = r
r+=0.05

endfor
for(m=0;m<=149;m+=1)

size1axis[m] = t
t+=1



endfor
duplicate FRET_Matrix_1Bin_Mass,ln_FRET_Matrix_1Bin_Mass
variable a,b
for(a=0;a<150;a+=1)
for(b=0;b<20;b+=1)
if(ln(FRET_Matrix_1Bin_Mass[a][b])>0)
ln_FRET_Matrix_1Bin_Mass[a][b]=ln(FRET_Matrix_1Bin_Mass[a][b])
else
ln_FRET_Matrix_1Bin_Mass[a][b]=0
endif
endfor
endfor

Display;AppendMatrixContour ln_FRET_Matrix_1Bin_Mass vs {FRETaxis,size1axis}
SetAxis bottom 0,1
SetAxis left 1,50
ModifyContour ln_FRET_Matrix_1Bin_Mass ctabLines={0,*,Rainbow256,1}
ModifyContour ln_FRET_Matrix_1Bin_Mass ctabLines={0,12,Rainbow256,1},autoLevels={0,12,100}
ModifyContour ln_FRET_Matrix_1Bin_Mass labels=0
ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
Label bottom "FRET Efficiency"
Label left "Size of Oligomer"
ModifyGraph width=226.772,height=226.772
ColorScale/C/N=text0/F=0/A=RC/E contour=ln_FRET_Matrix_1Bin_Mass

duplicate/o matrix_pascal matrix_pascal_mass

for(m=0;m<=(Dimsize(matrix_pascal,0));m+=1)

matrix_pascal_mass[m][0] = matrix_pascal_mass[m][0]
matrix_pascal_mass[m][1] = matrix_pascal_mass[m][1]*5
matrix_pascal_mass[m][2] = matrix_pascal_mass[m][2]*10
matrix_pascal_mass[m][3] = matrix_pascal_mass[m][3]*15
matrix_pascal_mass[m][4] = matrix_pascal_mass[m][4]*20
matrix_pascal_mass[m][5] = matrix_pascal_mass[m][5]*25
matrix_pascal_mass[m][6] = matrix_pascal_mass[m][6]*30
matrix_pascal_mass[m][7] = matrix_pascal_mass[m][7]*35
matrix_pascal_mass[m][8] = matrix_pascal_mass[m][8]*40
matrix_pascal_mass[m][9] = matrix_pascal_mass[m][9]*45
matrix_pascal_mass[m][10] = matrix_pascal_mass[m][10]*50
matrix_pascal_mass[m][11] = matrix_pascal_mass[m][11]*55
matrix_pascal_mass[m][12] = matrix_pascal_mass[m][12]*60
matrix_pascal_mass[m][13] = matrix_pascal_mass[m][13]*65
matrix_pascal_mass[m][14] = matrix_pascal_mass[m][14]*70
matrix_pascal_mass[m][15] = matrix_pascal_mass[m][15]*75
matrix_pascal_mass[m][16] = matrix_pascal_mass[m][16]*80
matrix_pascal_mass[m][17] = matrix_pascal_mass[m][17]*85
matrix_pascal_mass[m][18] = matrix_pascal_mass[m][18]*90
matrix_pascal_mass[m][19] = matrix_pascal_mass[m][19]*95
matrix_pascal_mass[m][20] = matrix_pascal_mass[m][20]*100
matrix_pascal_mass[m][21] = matrix_pascal_mass[m][21]*105
matrix_pascal_mass[m][22] = matrix_pascal_mass[m][22]*110
matrix_pascal_mass[m][23] = matrix_pascal_mass[m][23]*115
matrix_pascal_mass[m][24] = matrix_pascal_mass[m][24]*120
matrix_pascal_mass[m][25] = matrix_pascal_mass[m][25]*125
matrix_pascal_mass[m][26] = matrix_pascal_mass[m][26]*130
matrix_pascal_mass[m][27] = matrix_pascal_mass[m][27]*135
matrix_pascal_mass[m][28] = matrix_pascal_mass[m][28]*140
matrix_pascal_mass[m][29] = matrix_pascal_mass[m][29]*145

endfor
make/o/n=30 size5axis
variable u=0
for(m=0;m<=29;m+=1)
size5axis[m] = u
u+=5
endfor

killwaves matrix_pascal
duplicate/o matrix_pascal_mass FRET_Mass_Matrix
killwaves matrix_pascal_mass

//variable fact=A_bursts[301]
for(x=0;x<(dimsize(FRET_Matrix_1Bin,0));x+=1)
for(y=0;y<(dimsize(FRET_Matrix_1Bin,1));y+=1)
FRET_Matrix_1bin[x][y]=FRET_Matrix_1bin[x][y]
FRET_Matrix_1bin_mass[x][y]=FRET_Matrix_1bin_mass[x][y]
endfor
endfor
End 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to make histograms

function FRET_histograms()
wave variables,zmatrix,FRET_Matrix_1bin,FRETaxis,FRET_Small, FRET_Medium,FRET_Large,counter1,lnz_small,lnz_medium,lnz_large
variable n,m
make/o/n=20 FRET_Smallmers_Hist,FRET_Mediummers_Hist,FRET_Largemers_Hist,FRET_All
make/o/n=30 lnz_Smallmers_Hist,lnz_Mediummers_Hist,lnz_Largemers_Hist,lnz_All
print variables[2]
for(n=0;n<=variables[1];n+=1)
for(m=0;m<=19;m+=1)
FRET_Smallmers_Hist[m]=FRET_Smallmers_Hist[m]+FRET_Matrix_1bin[m][n]
endfor
endfor
for(n=(variables[1]+1);n<=variables[2];n+=1)
for(m=0;m<=19;m+=1)
FRET_Mediummers_Hist[m]=FRET_Mediummers_Hist[m]+FRET_Matrix_1bin[m][n]
endfor
endfor
for(n=(variables[2]+1);n<=150;n+=1)
for(m=0;m<=19;m+=1)
FRET_Largemers_Hist[m]=FRET_Largemers_Hist[m]+FRET_Matrix_1bin[m][n]
endfor
endfor



variable y = counter1[0]
for(m=0;m<=19;m+=1)

FRET_Small[m][y]=FRET_Smallmers_Hist[m]
FRET_Medium[m][y]=FRET_Mediummers_Hist[m]
FRET_Large[m][y]=FRET_Largemers_Hist[m]
FRET_ALL[m][y]=FRET_Smallmers_Hist[m]+FRET_Mediummers_Hist[m]+FRET_Largemers_Hist[m]
endfor



duplicate/o FRET_Small, root:FRET_Small
duplicate/o FRET_Medium, root:FRET_Medium
duplicate/o FRET_Large, root:FRET_Large
killwaves FRET_Small
killwaves FRET_Medium
killwaves FRET_Large

wave verylarge,data
variable t,k
make/o/n=1 larger=NAN,alloligomers=NAN
for(t=0;t<(dimsize(data,0));t+=1)
if(data[t][7]>150)
redimension/N=(k+1) larger
larger[k]=data[t][6]
k+=1
endif
endfor
make/o/n=20 verylarge_hist=NAN
Histogram/B={0,0.05,20} larger,verylarge_hist
for(m=0;m<=19;m+=1)

Verylarge[m][y]=verylarge_hist[m]
endfor
duplicate/o verylarge, root:verylarge
k=0
for(t=0;t<(dimsize(data,0));t+=1)
if(data[t][7]>0)
redimension/N=(k+1) alloligomers
alloligomers[k]=data[t][6]
k+=1
endif
endfor
wave allfret
make/o/n=20 alloligomers_hist=NAN
Histogram/B={0,0.05,20} alloligomers,alloligomers_hist
for(m=0;m<=19;m+=1)

allfret[m][y]=alloligomers_hist[m]
endfor
duplicate/o allfret, root:allfret
//Display FRET_Smallmers_Hist vs FRETaxis
//ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
//Label bottom "FRET Efficiency"
//Label left "#"
//TextBox/C/N=text0/F=0/A=MT/E "Smallmers"



//Display FRET_Mediummers_Hist vs FRETaxis
//ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
//Label bottom "FRET Efficiency"
//Label left "#"
//TextBox/C/N=text0/F=0/A=MT/E "Mediummers"
//Make/D/N=5/O W_coef_medium
//W_coef_medium[0] = {0.2,0.3,0.4,0.5,0.1,0.2}
//FuncFit/H="100100"/NTHR=0/TBOX=768 DoubleGauss W_coef_medium  FRET_mediummers_Hist /X=FRETaxis /D   

//Display FRET_Largemers_Hist vs FRETaxis
//ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
//Label bottom "FRET Efficiency"
//Label left "#"
//TextBox/C/N=text0/F=0/A=MT/E "Largemers"
//Make/D/N=6/O W_coef_large
//W_coef_large[0] = {0.2,0.3,0.4,0.4,0.1,0.2}
//FuncFit/H="100100"/NTHR=0/TBOX=768 DoubleGauss W_coef_large  FRET_largemers_Hist /X=FRETaxis /D   

End


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make ln(Ia/Id) matrices

function matrices_z()
wave data
variable l,count1,count2,k,count5,m
// Code to make 2D matrices
//Make the appropriate tables for adding data to
Make/O/N=(60,30) matrixz = 0
Make/O/N=(Dimsize(data,0),8) FretH

//Bin the size and FRET data
for(l=0;l<=(Dimsize(data,0));l+=1)
	FretH[l][0] = data[l][8]/0.1				// lnz
	FretH[l][1] = data[l][7]/5                 		// Size
	FretH[l][2] = round(FretH[l][0])			// Rounded FRET Efficiency
	FretH[l][3] = round(FretH[l][1])			// Rounded Size
	FretH[l][4] = FretH[l][0] - FretH[l][2]		// Needed to convert rounded to int
	FretH[l][5] = FretH[l][1] - FretH[l][3]		// Needed to convert rounded to int
if(FretH[l][4] < 0)
	FretH[l][6] = FretH[l][2] -1				// Ln(z) Int
else
	FretH[l][6] = FretH[l][2]				// Ln(z) Efficiency Int
endif
if(FretH[l][5] < 0)
	FretH[l][7] = FretH[l][3] -1				// Size Int
else
	FretH[l][7] = FretH[l][3]				// Size Int
endif
endfor

////The following is data not accounting for Pascal

//Delete any data in which the size doesn't fit into the matrix
for(k=0;k<=(Dimsize(FretH,0));k+=1)
		if(FretH[k][6] > 30)
		DeletePoints k,1, FretH
		count1 +=1
		endif
endfor
for(k=0;k<=(Dimsize(FretH,0));k+=1)
		if(FretH[k][6] < -30)
		DeletePoints k,1, FretH
		count1 +=1
		endif
endfor
for(k=0;k<=(Dimsize(FretH,0));k+=1)
	if(FretH[k][7] > 30)
	DeletePoints k,1, FretH
	count2 +=1
	endif
endfor
Printf "%d data points deleted due to FRET being too high.\r",count1
Printf "%d data points deleted due to Size being too high.\r",count2

//Add the data to the matrixz.
variable x,y
for(m=0;m<=(Dimsize(FretH,0));m+=1)
if(FretH[m][7]<=29)
if(FretH[m][6]<=29)
if(FretH[m][6]>-29)
y = FretH[m][6]+30
x = FretH[m][7] 
duplicate/o matrixz, matrixz10
matrixz[y][x]=matrixz10[y][x]+1 // Adds 1 each time data needs to be added
endif
endif
endif
endfor

////Accounting for Pascal
Make/O/N=(60,150) pascalm = 0
Make/O/N=(Dimsize(data,0),8) Pascal_Fret
for(l=0;l<=(Dimsize(data,0));l+=1)
	Pascal_Fret[l][0] = data[l][8]/0.1				// Bin lnz
	Pascal_Fret[l][1] = data[l][7]					//sizes not binned. 
	Pascal_Fret[l][2] = round(Pascal_Fret[l][0])
	Pascal_Fret[l][3] = round(Pascal_Fret[l][1])
	Pascal_Fret[l][4] = Pascal_Fret[l][0] - Pascal_Fret[l][2]
	Pascal_Fret[l][5] = Pascal_Fret[l][1] - Pascal_Fret[l][3]
if(Pascal_Fret[l][4] < 0)
	Pascal_Fret[l][6] = Pascal_Fret[l][2] -1
else
	Pascal_Fret[l][6] = Pascal_Fret[l][2]		
endif
if(Pascal_FRET[l][5] < 0)
	Pascal_Fret[l][7] = Pascal_Fret[l][3] -1
else
	Pascal_Fret[l][7] = Pascal_Fret[l][3]		// Sizes
endif
endfor

//Delete points that won't fit into the 150x60 matrixz

for(k=0;k<=(Dimsize(Pascal_Fret,0));k+=1)
		if(Pascal_Fret[k][6] > 30)
		DeletePoints k,1, Pascal_Fret
		endif
endfor
for(k=0;k<=(Dimsize(Pascal_Fret,0));k+=1)
		if(Pascal_Fret[k][6] <-30)
		DeletePoints k,1, Pascal_Fret
		endif
endfor

for(k=0;k<=(Dimsize(Pascal_Fret,0));k+=1)
	if(FretH[k][7] > 150)
	DeletePoints k,1, Pascal_Fret
	endif

endfor

//Put data into the 150x60 matrixz

for(m=0;m<=(Dimsize(Pascal_Fret,0));m+=1)
variable count4
if( Pascal_Fret[m][7]<=149)
if(Pascal_Fret[m][6]<=29)
if(Pascal_Fret[m][6]>-29)
y = Pascal_Fret[m][6]+30
x = Pascal_Fret[m][7] 
duplicate/o pascalm, pascalm10
pascalm[y][x]=pascalm10[y][x]+1
endif
endif
endif

endfor


//Multiply size distributions by Pascal factor to account for "invisible" oligomers

for(m=0;m<=(Dimsize(pascalm,0));m+=1)

//pascalm[m][2]=pascalm[m][2]*2
//pascalm[m][3]=pascalm[m][3]*1.33
//pascalm[m][4]=pascalm[m][4]*1.14
//pascalm[m][5]=pascalm[m][5]*1.067

endfor


////Now to modify the matrixz not accounted for Pascal
duplicate/o matrixz, matrixz_pascal
for(m=0;m<=(Dimsize(matrixz_pascal,0));m+=1)

matrixz_pascal[m][0]=pascalm[m][0]+pascalm[m][1]+pascalm[m][2]+pascalm[m][3]+pascalm[m][4]
matrixz_pascal[m][1]=pascalm[m][5]+pascalm[m][6]+pascalm[m][7]+pascalm[m][8]+pascalm[m][9]

endfor
duplicate/o pascalm lnz_Matrix_1bin
killwaves pascalm
killwaves pascalm10
//killwaves Pascal_FRET
killwaves matrixz10
killwaves FretH
killwaves matrixz
duplicate matrixz_pascal lnz_Matrix

duplicate/o matrixz_pascal matrixz_pascal_mass
duplicate/o lnz_Matrix_1Bin lnz_Matrix_1Bin_Mass
variable n=1
for(k=0;k<=(Dimsize(lnz_Matrix_1Bin,1));k+=1)
for(m=0;m<=(Dimsize(lnz_Matrix_1Bin,0));m+=1)

lnz_Matrix_1Bin_Mass[m][k] = lnz_Matrix_1Bin[m][k]*k


endfor
endfor
make/o/n=60 lnzaxis
make/o/n=150 sizeaxis

variable r = -2.95
variable t = 0
variable u=0
for(m=0;m<=59;m+=1)

lnzaxis[m] = r
r+=0.1

endfor
for(m=0;m<=149;m+=1)

sizeaxis[m] = u
u+=1

endfor



Display;AppendMatrixContour lnz_Matrix_1Bin_Mass vs {lnzaxis,sizeaxis}
SetAxis bottom -3,3
SetAxis left 1,50
ModifyContour lnz_Matrix_1Bin_Mass labels=0
ModifyContour lnz_Matrix_1Bin_Mass autoLevels={*,*,50}
ModifyContour lnz_Matrix_1Bin_Mass ctabLines={0,*,Rainbow256,1}
ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
Label bottom "ln(Ia/Id)"
Label left "Size of Oligomer"
ModifyGraph width=226.772,height=226.772
ColorScale/C/N=text0/F=0/A=RC/E contour=lnz_Matrix_1Bin_Mass
for(m=0;m<=(Dimsize(matrixz_pascal,0));m+=1)

matrixz_pascal_mass[m][0] = matrixz_pascal_mass[m][0]
matrixz_pascal_mass[m][1] = matrixz_pascal_mass[m][1]*5
matrixz_pascal_mass[m][2] = matrixz_pascal_mass[m][2]*10
matrixz_pascal_mass[m][3] = matrixz_pascal_mass[m][3]*15
matrixz_pascal_mass[m][4] = matrixz_pascal_mass[m][4]*20
matrixz_pascal_mass[m][5] = matrixz_pascal_mass[m][5]*25
matrixz_pascal_mass[m][6] = matrixz_pascal_mass[m][6]*30
matrixz_pascal_mass[m][7] = matrixz_pascal_mass[m][7]*35
matrixz_pascal_mass[m][8] = matrixz_pascal_mass[m][8]*40
matrixz_pascal_mass[m][9] = matrixz_pascal_mass[m][9]*45
matrixz_pascal_mass[m][10] = matrixz_pascal_mass[m][10]*50
matrixz_pascal_mass[m][11] = matrixz_pascal_mass[m][11]*55
matrixz_pascal_mass[m][12] = matrixz_pascal_mass[m][12]*60
matrixz_pascal_mass[m][13] = matrixz_pascal_mass[m][13]*65
matrixz_pascal_mass[m][14] = matrixz_pascal_mass[m][14]*70
matrixz_pascal_mass[m][15] = matrixz_pascal_mass[m][15]*75
matrixz_pascal_mass[m][16] = matrixz_pascal_mass[m][16]*80
matrixz_pascal_mass[m][17] = matrixz_pascal_mass[m][17]*85
matrixz_pascal_mass[m][18] = matrixz_pascal_mass[m][18]*90
matrixz_pascal_mass[m][19] = matrixz_pascal_mass[m][19]*95
matrixz_pascal_mass[m][20] = matrixz_pascal_mass[m][20]*100
matrixz_pascal_mass[m][21] = matrixz_pascal_mass[m][21]*105
matrixz_pascal_mass[m][22] = matrixz_pascal_mass[m][22]*110
matrixz_pascal_mass[m][23] = matrixz_pascal_mass[m][23]*115
matrixz_pascal_mass[m][24] = matrixz_pascal_mass[m][24]*120
matrixz_pascal_mass[m][25] = matrixz_pascal_mass[m][25]*125
matrixz_pascal_mass[m][26] = matrixz_pascal_mass[m][26]*130
matrixz_pascal_mass[m][27] = matrixz_pascal_mass[m][27]*135
matrixz_pascal_mass[m][28] = matrixz_pascal_mass[m][28]*140
matrixz_pascal_mass[m][29] = matrixz_pascal_mass[m][29]*145

endfor
killwaves matrixz_pascal
duplicate/o matrixz_pascal_mass lnz_Mass_Matrix
killwaves matrixz_pascal_mass
End

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to make lnz histograms for each size boundary. 
function lnz_histograms()
wave variables, lnz_Matrix_1bin,lnzaxis,lnz_Small, lnz_Medium,lnz_Large
variable n,m
make/o/n=60 lnz_Smallmers_Hist,lnz_Mediummers_Hist,lnz_Largemers_Hist
for(n=0;n<=variables[1];n+=1)
for(m=0;m<=59;m+=1)
lnz_Smallmers_Hist[m]=lnz_Smallmers_Hist[m]+lnz_Matrix_1bin[m][n]
endfor
endfor
for(n=(variables[1]+1);n<=variables[2];n+=1)
for(m=0;m<=59;m+=1)
lnz_Mediummers_Hist[m]=lnz_Mediummers_Hist[m]+lnz_Matrix_1bin[m][n]
endfor
endfor
for(n=(variables[2]+1);n<=150;n+=1)
for(m=0;m<=59;m+=1)
lnz_Largemers_Hist[m]=lnz_Largemers_Hist[m]+lnz_Matrix_1bin[m][n]
endfor
endfor
wave counter1
variable y = counter1[0]
for(m=0;m<=59;m+=1)

lnz_Small[m][y]=lnz_Smallmers_Hist[m]
lnz_Medium[m][y]=lnz_Mediummers_Hist[m]
lnz_Large[m][y]=lnz_Largemers_Hist[m]
endfor
duplicate/o lnz_Small, root:lnz_Small
duplicate/o lnz_Medium, root:lnz_Medium
duplicate/o lnz_Large, root:lnz_Large
killwaves lnz_Small
killwaves lnz_Medium
killwaves lnz_Large

//Display lnz_Smallmers_Hist vs lnzaxis
//ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
//Label bottom "lnz Efficiency"
//Label left "#"
//TextBox/C/N=text0/F=0/A=MT/E "Smallmers"

//Display lnz_Mediummers_Hist vs lnzaxis
//ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
//Label bottom "lnz Efficiency"
//Label left "#"
//TextBox/C/N=text0/F=0/A=MT/E "Mediummers"

//Display lnz_Largemers_Hist vs lnzaxis
//ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
//Label bottom "lnz Efficiency"
//Label left "#"
//TextBox/C/N=text0/F=0/A=MT/E "Largemers"




End


function zplot()
wave data
make/o/n=(30,150) zmatrix
make/o/n=30 tempz

variable a,b,c,d,e

for(a=0;a<7;a+=1)
c=0
make/o/n=1 temp
for(b=0;b<(dimsize(data,0));b+=1)
if(data[b][7]>a && data[b][7]<(a+2))
redimension/n=(c) temp
temp[c]=data[b][8]
c+=1
endif

endfor
Make/N=30/O temp_Hist,z_axis;DelayUpdate
Histogram/B={-3,0.2,30} temp,temp_Hist
e=-3
for(d=0;d<30;d+=1)
zmatrix[d][a]=temp_hist[d]
z_axis[d]=e
e+=0.2
endfor
killwaves temp,temp_hist
endfor
end

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to generate printable layout of data. 

Function see()

wave FRET_matrix_1bin,FRETaxis,size1axis,hour,A_bursts

duplicate/o FRET_Matrix_1Bin,ln_FRET_Matrix_1Bin
variable a,b
for(a=0;a<150;a+=1)
for(b=0;b<20;b+=1)
if(ln(FRET_Matrix_1Bin[a][b])>0)
ln_FRET_Matrix_1Bin[a][b]=ln(FRET_Matrix_1Bin[a][b])
else
ln_FRET_Matrix_1Bin[a][b]=0
endif
endfor
endfor


duplicate/o FRET_Matrix_1Bin,see_lnFRETMatrix

for(a=0;a<150;a+=1)
for(b=0;b<20;b+=1)
if(ln(FRET_Matrix_1Bin[a][b])>0)
see_lnFRETMatrix[a][b]=ln(FRET_Matrix_1Bin[a][b])
else
see_lnFRETMatrix[a][b]=0
endif
endfor
endfor

Display;AppendMatrixContour see_lnFRETMatrix vs {FRETaxis,size1axis}
SetAxis bottom 0,1
SetAxis left 1,50
ModifyContour see_lnFRETMatrix labels=0
ModifyContour see_lnFRETMatrix ctabLines={0,12,Rainbow256,1},autoLevels={0,12,50}
ColorScale/C/N=text0 "ln(counts)"
ModifyContour see_lnFRETMatrix autoLevels={*,*,50}
ModifyContour see_lnFRETMatrix ctabLines={0,*,Rainbow256,1}
ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
Label bottom "FRET Efficiency"
Label left "Size of Oligomer"
ModifyGraph width=100,height=100
ColorScale/C/N=text0/F=0/A=RC/E contour=see_lnFRETMatrix
TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=0.00/E "Number Distribution"
ModifyGraph swapXY=1
ModifyGraph manTick(bottom)={0,10,0,0},manMinor(bottom)={0,0}
ModifyContour see_lnFRETMatrix ctabLines={0,7,Rainbow256,1},autoLevels={0,7,50}

//wave lnz_Matrix_1Bin_Mass,lnzaxis,sizeaxis
//Display;AppendMatrixContour lnz_Matrix_1Bin_Mass vs {lnzaxis,sizeaxis}
//SetAxis bottom -3,3
//SetAxis left 1,50
//ModifyContour lnz_Matrix_1Bin_Mass labels=0
//ModifyContour lnz_Matrix_1Bin_Mass autoLevels={*,*,50}
//ModifyContour lnz_Matrix_1Bin_Mass ctabLines={0,*,Rainbow256,1}
//ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
//Label bottom "ln(Ia/Id)"
//Label left "Size of Oligomer"
//ModifyGraph width=226.772,height=226.772
//ColorScale/C/N=text0/F=0/A=RC/E contour=lnz_Matrix_1Bin_Mass


//wave lnz_smallmers, lnz_mediummers,lnz_largemers

//Display lnz_Smallmers_Hist vs lnzaxis
//ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
//Label bottom "lnz Efficiency"
//Label left "#"
//TextBox/C/N=text0/F=0/A=MT/E "Smallmers"

//Display lnz_Mediummers_Hist vs lnzaxis
//ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
//Label bottom "lnz Efficiency"
//Label left "#"
//TextBox/C/N=text0/F=0/A=MT/E "Mediummers"

//Display lnz_Largemers_Hist vs lnzaxis
//ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
//Label bottom "lnz Efficiency"
//Label left "#"
//TextBox/C/N=text0/F=0/A=MT/E "Largemers"

wave FRET_smallmers_hist, FRET_mediummers_hist,FRET_largemers_hist,FRET_ALL,donor_bursts

make/o/n=1 Small_Tot, Medium_Tot, Large_tot
variable n
for(n=0;n<=19;n+=1)
small_tot[0]=small_tot[0]+Fret_smallmers_hist[n]
medium_tot[0]=medium_tot[0]+Fret_mediummers_hist[n]
large_tot[0]=large_tot[0]+Fret_largemers_hist[n]
endfor
variable h=hour[0]
wave A_rate
string h1=num2str(h)
variable s = round((small_tot[0]))
variable m = round((medium_tot[0]))
variable l = round((large_tot[0]))
variable tot=s+m+l
variable mon=donor_bursts[0]
variable sp = s
variable mp = m
variable lp = l
string s1=num2str(s)
string m1=num2str(m)
string l1=num2str(l)
string tot1=num2str(tot)
string sp1=num2str(sp)
string mp1=num2str(mp)
string lp1=num2str(lp)
string monom=num2str(mon)

Display FRET_Smallmers_Hist vs FRETaxis
ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
Label bottom "FRET Efficiency"
Label left "# events"
TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=0.00/E h1+" hours Smallmers"
TextBox/C/N=text0/F=0/A=MB/X=0.00/Y=0.00/E "\Z08\K(65280,0,0)"+s1+" out of "+tot1+" events"

Display FRET_Mediummers_Hist vs FRETaxis
ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
Label bottom "FRET Efficiency"
Label left "# events"
TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=0.00/E h1+" hours Mediummers"
TextBox/C/N=text0/F=0/A=MB/X=0.00/Y=0.00/E "\Z08\K(65280,0,0)"+m1+" out of "+tot1+" events"

Display FRET_Largemers_Hist vs FRETaxis
ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
Label bottom "FRET Efficiency"
Label left "# events"
TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=0.00/E h1+"  hours Largemers"
TextBox/C/N=text0/F=0/A=MB/X=0.00/Y=0.00/E "\Z08\K(65280,0,0)"+l1+" out of "+tot1+" events"
wave stats
display stats[][0]
TextBox/C/N=text0/F=0/A=MT/E Pathway
Label bottom "Frame Number"
Label left "Event rate (s\\S-1\\M)"
AppendToGraph/R stats[][1]
Label right "Event rate (s\\S-1\\M)"
ModifyGraph rgb(Stats#1)=(0,0,65280)
ModifyGraph axRGB(right)=(0,0,65280)
ModifyGraph axRGB(left)=(65280,3328,0),tlblRGB(left)=(65280,3328,0);DelayUpdate
ModifyGraph tlblRGB(right)=(0,0,65280),alblRGB(left)=(65280,3328,0);DelayUpdate
ModifyGraph alblRGB(right)=(0,0,65280)
SetAxis left 0,*
SetAxis right 0,*
Legend/C/N=text1/J/A=RC/E "\\s(Stats) Donor Bursts\r\\s(Stats#1) Coincident Bursts"

// wave Division,A_rate,B_rate,C_rate
//Display A_rate vs Division
//Label bottom "Frame Number"
//Label left "Burstrate"
//TextBox/C/N=text0/F=0/A=MT/X=0.00/Y=0.00/E "Donor Burstrate"
//ModifyGraph manTick(bottom)={0,100,0,0},manMinor(bottom)={0,0}

//Display B_rate vs Division
//Label bottom "Frame Number"
//Label left "Burstrate"
//TextBox/C/N=text0/F=0/A=MT/X=0.00/Y=0.00/E "Acceptor Burstrate"
//ModifyGraph manTick(bottom)={0,100,0,0},manMinor(bottom)={0,0}

//Display C_rate vs Division
//Label bottom "Frame Number"
//Label left "Burstrate"
//TextBox/C/N=text0/F=0/A=MT/exX=0.00/Y=0.00/E "Coincident Burstrate"
//ModifyGraph manTick(bottom)={0,100,0,0},manMinor(bottom)={0,0}
wave toobig
string toobig1=num2str(toobig[0])
wave mb1
string mb2=num2str(mb1[0])

newLayout/N=p
textBox/C/N=text1/F=0/A=MT/X=0.00/Y=0/E h1+" Hours "
AppendLayoutObject/r=(72,90,302,240)/F=1 graph graph0
AppendLayoutObject/r=(72,260,232,370)/F=1 graph graph1
AppendLayoutObject/r=(221,260,372,370)/F=1 graph graph2
AppendLayoutObject/r=(371,260,522,370)/F=1 graph graph3

drawline 70,385,522,385
DrawRect 328.5,117,486.75,131.25
DrawText 391.5,129.75,"\Z08 Statistics"
DrawRect 328.5,131.25,407.625,145.25
DrawRect 407.625,131.25,486.75,145.25
DrawText 333,144,"\\Z08Oligomer Size > 150"
DrawText 417,144,"\\Z08"+toobig1+" Oligomers"
DrawRect 328.5,145.25,407.625,159.25
DrawRect 407.625,145.25,486.75,159.25
DrawText 331.5,157.75,"\\Z08 Useful detected"
DrawText 417,157.75,"\\Z08"+tot1+"/"+monom+""
DrawRect 328.5,159.25,407.625,173.25
DrawRect 407.625,159.25,486.75,173.25
DrawText 331.5,171.5,"\\Z08 Monomer brightness"
DrawText 417,171.5,"\\Z08"+mb2

AppendLayoutObject/r=(72,400,522,510)/F=1 graph graph4
//AppendLayoutObject/r=(221,400,372,510)/F=1 graph graph5
//AppendLayoutObject/r=(371,400,522,510)/F=1 graph graph6
//SetDrawEnv dash= 1,fillpat= 0;DelayUpdate
//DrawRect 72,72,521,522

end

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to generate FRET histogram for each size. 

function singlehist(path)
variable path
string p
p=num2str(path)
setdatafolder root:$p
wave FRET_Matrix_1Bin,FRETaxis
variable n,m
for(n=0;n<=20;n+=1)
Display FRET_Matrix_1Bin[][n] vs FRETaxis 
ModifyGraph mode=5,hbFill=3,rgb=(0,0,0)
Label bottom "FRET Efficiency"
Label left "#"
string r = num2str(n)
TextBox/C/N=text0/F=0/A=MT/E "Oligomer Size = "+r 
endfor

newLayout/N=Histograms_1
AppendLayoutObject/r=(75,95,280,230)/F=1 graph graph0
AppendLayoutObject/r=(75,230,280,365)/F=1 graph graph2
AppendLayoutObject/r=(75,365,280,500)/F=1 graph graph4
AppendLayoutObject/r=(75,500,280,635)/F=1 graph graph6
AppendLayoutObject/r=(75,635,280,770)/F=1 graph graph8
AppendLayoutObject/r=(315,95,520,230)/F=1 graph graph1
AppendLayoutObject/r=(315,230,520,365)/F=1 graph graph3
AppendLayoutObject/r=(315,365,520,500)/F=1 graph graph5
AppendLayoutObject/r=(315,500,520,635)/F=1 graph graph7
AppendLayoutObject/r=(315,635,520,770)/F=1 graph graph9
TextBox/C/N=text0/A=MT/X=0.00/Y=0/f=0 "\\f01\\Z12"+p+"\\f01\\Z16 Hours"
newLayout/N=Histograms_2
AppendLayoutObject/r=(75,95,280,230)/F=1 graph graph10
AppendLayoutObject/r=(75,230,280,365)/F=1 graph graph12
AppendLayoutObject/r=(75,365,280,500)/F=1 graph graph14
AppendLayoutObject/r=(75,500,280,635)/F=1 graph graph16
AppendLayoutObject/r=(75,635,280,770)/F=1 graph graph18
AppendLayoutObject/r=(315,95,520,230)/F=1 graph graph11
AppendLayoutObject/r=(315,230,520,365)/F=1 graph graph13
AppendLayoutObject/r=(315,365,520,500)/F=1 graph graph15
AppendLayoutObject/r=(315,500,520,635)/F=1 graph graph17
AppendLayoutObject/r=(315,635,520,770)/F=1 graph graph19
TextBox/C/N=text0/A=MT/X=0.00/Y=0/f=0 "\\f01\\Z12"+p+"\\f01\\Z16 Hours"
//AppendLayoutObject graph graph2
end



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Layout all 2D number distributions on one page. 

function print2d(z)
variable z
variable n,m
for(n=1;n<=z;n+=1)
string p
p=num2str(n)
setdatafolder root:$p
wave ln_FRET_matrix_1bin_mass,FRETaxis,size1axis
Display;AppendMatrixContour ln_FRET_Matrix_1Bin_Mass vs {FRETaxis,size1axis}
SetAxis bottom 0,1
SetAxis left 1,50
ModifyContour ln_FRET_Matrix_1Bin_Mass labels=0
ModifyContour ln_FRET_Matrix_1Bin_Mass autoLevels={*,*,50}
ModifyContour ln_FRET_Matrix_1Bin_Mass ctabLines={0,12,Rainbow256,1},autoLevels={0,12,50}
ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
Label bottom "FRET Efficiency"
Label left "Size of Oligomer"
colorScale/C/N=text1/F=0/A=RC/E contour=ln_FRET_Matrix_1Bin_Mass
TextBox/C/N=text0/F=0/A=MT/E "Timepoint "+p
wave lnz_Matrix_1Bin_Mass,lnzaxis,sizeaxis
Display;AppendMatrixContour lnz_Matrix_1Bin_Mass vs {lnzaxis,sizeaxis}
SetAxis bottom -3,3
SetAxis left 1,50
ModifyContour lnz_Matrix_1Bin_Mass labels=0
ModifyContour lnz_Matrix_1Bin_Mass autoLevels={*,*,50}
ModifyContour lnz_Matrix_1Bin_Mass ctabLines={0,*,Rainbow256,1}
ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
Label bottom "ln(Ia/Id)"
Label left "Size of Oligomer"
colorScale/C/N=text1/F=0/A=RC/E contour=lnz_Matrix_1Bin_Mass
TextBox/C/N=text0/F=0/A=MT/E "Timepoint "+p
endfor
newLayout/N=Two_Dimensional_1
AppendLayoutObject/r=(75,85,280,220)/F=1 graph graph0
AppendLayoutObject/r=(75,220,280,355)/F=1 graph graph2
AppendLayoutObject/r=(75,365,280,490)/F=1 graph graph4
AppendLayoutObject/r=(75,490,280,625)/F=1 graph graph6
AppendLayoutObject/r=(75,625,280,760)/F=1 graph graph8
AppendLayoutObject/r=(315,85,520,220)/F=1 graph graph1
AppendLayoutObject/r=(315,220,520,355)/F=1 graph graph3
AppendLayoutObject/r=(315,355,520,490)/F=1 graph graph5
AppendLayoutObject/r=(315,490,520,625)/F=1 graph graph7
AppendLayoutObject/r=(315,625,520,760)/F=1 graph graph9
newLayout/N=Two_Dimensional_2
AppendLayoutObject/r=(75,75,280,210)/F=1 graph graph10
AppendLayoutObject/r=(75,210,280,345)/F=1 graph graph12
AppendLayoutObject/r=(75,345,280,480)/F=1 graph graph14
AppendLayoutObject/r=(75,480,280,615)/F=1 graph graph16
AppendLayoutObject/r=(75,615,280,750)/F=1 graph graph18
AppendLayoutObject/r=(315,75,520,210)/F=1 graph graph11
AppendLayoutObject/r=(315,210,520,345)/F=1 graph graph13
AppendLayoutObject/r=(315,345,520,480)/F=1 graph graph15
AppendLayoutObject/r=(315,480,520,615)/F=1 graph graph17
AppendLayoutObject/r=(315,615,520,750)/F=1 graph graph19


//AppendLayoutObject graph graph2
end

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Totals over time. 

function totalsFRET()
setdatafolder root:
variable n
variable m
variable smalltot,mediumlowtot,mediumhightot,largelowtot,largehightot
wave timepoints,Fret_Large,Fret_Medium,Fret_small
make/o/n=(dimsize(timepoints,0)) Small,largelow,largehigh, MediumLow, MediumHigh, Large,all,SmallNorm,MediumLowNorm, MediumHighNorm,LargeNorm,Monomerb,speca,specb
for(n=0;n<=(dimsize(timepoints,0));n+=1)
smalltot=0
mediumlowtot=0
mediumhightot=0
largelowtot=0
largehightot=0
for(m=0;m<=19;m+=1)
smalltot=smalltot+Fret_small[m][n+1]

endfor
for(m=0;m<=7;m+=1)
mediumlowtot=mediumlowtot+Fret_medium[m][n+1]
largelowtot=largelowtot+Fret_large[m][n+1]
endfor
for(m=8;m<=19;m+=1)
mediumhightot=mediumhightot+Fret_medium[m][n+1]
largehightot=largehightot+Fret_large[m][n+1]
endfor

small[n]=smalltot
mediumlow[n]=mediumlowtot
mediumhigh[n]=mediumhightot
largehigh[n]=largehightot
largelow[n]=largelowtot

all[n]=(smalltot+mediumlowtot+largelowtot+largehightot+mediumhightot)
speca[n]=(smalltot+mediumlowtot+largelowtot)
specb[n]=(mediumhightot+largehightot)
endfor
display small vs timepoints
AppendToGraph Large vs timepoints
AppendToGraph Mediumlow  vs timepoints
AppendToGraph Mediumhigh  vs timepoints
AppendToGraph largelow  vs timepoints
AppendToGraph largehigh  vs timepoints
AppendToGraph All vs timepoints
Label bottom "Time / hours"
ModifyGraph rgb(Small)=(65280,43520,0),rgb(Largelow)=(0,0,65280);DelayUpdate
ModifyGraph rgb(largehigh)=(39168,13056,0)
ModifyGraph rgb(Mediumlow)=(0,52224,0)
ModifyGraph rgb(MediumHigh)=(29440,0,58880)
Legend/C/N=text0/J/F=0/A=MT/E "\\s(Small) Smallmers\r\\s(MediumLow) Mediummers (low FRET)\r\\s(MediumHigh) Mediummers (high FRET)\r\\s(Largelow) Largemers (Low FRET)\r\\s(Largehigh) Largemers (High FRET)\r\\s(all) All"
Legend/C/N=text0/J/A=RC/X=0.00/Y=0.00
ModifyGraph mode(all)=3,marker(all)=19
ModifyGraph msize(all)=1.5
ModifyGraph mode(Small)=3,marker(Small)=19,msize(Small)=1.5
ModifyGraph mode(MediumHigh)=4,marker(MediumHigh)=19,msize(MediumHigh)=1.5;DelayUpdate
ModifyGraph mode(all)=4
ModifyGraph mode(Small)=4
ModifyGraph mode(MediumLow)=4,marker(MediumLow)=19,msize(MediumLow)=1.5
ModifyGraph mode=4,marker=19,msize=1.5
Label left "Number of Molecules"

display speca,specb vs timepoints
Label bottom "Time / hours"
Label left "Number of Molecules"
ModifyGraph mode=4,marker=19,msize=1.5,rgb(speca)=(0,0,65280)
Legend/C/N=text1/F=0/A=RC/E "\\s(speca) Low FRET\r\\s(specb) High FRET"

//display smallnorm vs timepoints
//AppendToGraph Largenorm vs timepoints
//AppendToGraph Mediumlownorm  vs timepoints
//AppendToGraph Mediumhighnorm  vs timepoints
//Label bottom "Time / hours"
//Label left "Fraction of Oligomers"
//ModifyGraph rgb(Smallnorm)=(65280,43520,0),rgb(Largenorm)=(0,0,65280);DelayUpdate
//ModifyGraph rgb(Mediumlownorm)=(0,52224,0)
//ModifyGraph rgb(Mediumhighnorm)=(29440,0,58880)
//Legend/C/N=text0/J/F=0/A=MT/E "\\s(Smallnorm) Smallmers\r\\s(Largenorm) Largemers\r\\s(MediumLownorm) Mediummers (low FRET)\r\\s(MediumHighnorm) Mediummers (high FRET)"
//Legend/C/N=text0/J/A=RC/X=0.00/Y=0.00


end

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// All 2D layouts. 
function copy_plots()
wave  FRET_Matrix_1Bin,fretaxis,size1axis

duplicate/o FRET_Matrix_1Bin,seelnFRET_Matrix_1Bin
variable a,b
for(a=0;a<150;a+=1)
for(b=0;b<20;b+=1)
if(ln(FRET_Matrix_1Bin[a][b])>0)
seelnFRET_Matrix_1Bin[a][b]=ln(FRET_Matrix_1Bin[a][b])
else
seelnFRET_Matrix_1Bin[a][b]=0
endif
endfor
endfor

end
function twod()
kill()
setdatafolder root:
wave timepoints
variable n
for(n=0;n<(dimsize(timepoints,0));n+=1)
setdatafolder root:
variable l=timepoints[n]
string y=num2str(l)
setdatafolder $y
copy_plots()

wave seelnFRET_Matrix_1Bin


// scale etc
Display;AppendMatrixContour seelnFRET_Matrix_1Bin vs {FRETaxis,size1axis}
SetAxis bottom 0,1
SetAxis left 1,150
ModifyContour seelnFRET_Matrix_1Bin labels=0
ModifyContour seelnFRET_Matrix_1Bin autoLevels={*,*,50}
ModifyContour seelnFRET_Matrix_1Bin ctabLines={0,12,Rainbow256,1},autoLevels={0,12,50}
ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
ModifyGraph manTick(left)={0,20,0,0},manMinor(left)={0,0}
ModifyGraph fSize(bottom)=8
ModifyGraph fSize(left)=8
ModifyGraph swapXY=1
ModifyGraph width=120
ModifyGraph height=100
ModifyGraph noLabel=2
TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=0.00/E "\Z08"+y+" hours number distribution"
TextBox/C/N=text1/X=8.00/Y=0.00
ModifyContour seelnFRET_Matrix_1Bin autoLevels={0,5,50}
SetAxis bottom 1,50
ModifyGraph btLen=5,manTick(bottom)={0,10,0,0},manMinor(bottom)={0,0}
ModifyGraph mirror=1
ModifyContour seelnFRET_Matrix_1Bin ctabLines={0,5,Rainbow256,1}
endfor
Display;AppendMatrixContour seelnFRET_Matrix_1Bin vs {FRETaxis,size1axis}
SetAxis bottom 0,1
SetAxis left 1,150
ModifyContour seelnFRET_Matrix_1Bin labels=0
ModifyContour seelnFRET_Matrix_1Bin autoLevels={*,*,50}
ModifyContour seelnFRET_Matrix_1Bin ctabLines={0,12,Rainbow256,1},autoLevels={0,12,50}
ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
Label bottom "\Z08 FRET Efficiency"
Label left "\Z08 Size of Oligomer"
ColorScale/C/N=text4/F=0/A=RC/E frame=0.00
TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=0.00/E "Scale bar etc"
ModifyGraph manTick(left)={0,20,0,0},manMinor(left)={0,0}
ModifyGraph fSize(bottom)=8
ModifyGraph fSize(left)=8
ModifyGraph swapXY=1
ModifyContour seelnFRET_Matrix_1Bin autoLevels={0,5,50}
SetAxis bottom 1,50
ModifyGraph btLen=5,manTick(bottom)={0,10,0,0},manMinor(bottom)={0,0}
ModifyGraph mirror=1
ModifyContour seelnFRET_Matrix_1Bin ctabLines={0,5,Rainbow256,1}

newLayout/N=Histograms_1
AppendLayoutObject/r=(75,95,280,230)/F=0 graph graph0
AppendLayoutObject/r=(75,229,280,364)/F=0 graph graph2
AppendLayoutObject/r=(75,363,280,498)/F=0 graph graph4
AppendLayoutObject/r=(75,497,280,632)/F=0 graph graph6
AppendLayoutObject/r=(75,631,280,766)/F=0 graph graph8
AppendLayoutObject/r=(279,95,484,230)/F=0 graph graph1
AppendLayoutObject/r=(279,229,484,364)/F=0 graph graph3
AppendLayoutObject/r=(279,363,484,498)/F=0 graph graph5
AppendLayoutObject/r=(279,497,484,632)/F=0 graph graph7
AppendLayoutObject/r=(279,631,484,766)/F=0 graph graph9

newLayout/N=Histograms_2
AppendLayoutObject/r=(75,95,280,230)/F=0 graph graph10
AppendLayoutObject/r=(75,229,280,364)/F=0 graph graph12
AppendLayoutObject/r=(75,363,280,498)/F=0 graph graph14
AppendLayoutObject/r=(75,497,280,632)/F=0 graph graph16
AppendLayoutObject/r=(75,631,280,766)/F=0 graph graph18
AppendLayoutObject/r=(279,95,484,230)/F=0 graph graph11
AppendLayoutObject/r=(279,229,484,364)/F=0 graph graph13
AppendLayoutObject/r=(279,363,484,498)/F=0 graph graph15
AppendLayoutObject/r=(279,497,484,632)/F=0 graph graph17
AppendLayoutObject/r=(279,631,484,766)/F=0 graph graph19


//AppendLayoutObject graph graph2
end


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Layout mass plots.


function twodmass()
setdatafolder root:
wave timepoints
variable n
for(n=0;n<(dimsize(timepoints,0));n+=1)
setdatafolder root:
variable l=timepoints[n]
string y=num2str(l)
setdatafolder $y
wave  FRET_Matrix_1bin_mass,fretaxis,size1axis
make/o/n=(20,150) FRET_Matrix_1bin_mass_conc
variable s,t
for(s=0;s<=19;s+=1)
for(t=0;t<=149;t+=1)
FRET_Matrix_1bin_mass_conc[s][t]=70*FRET_Matrix_1bin_mass[s][t]
endfor
endfor
Display;AppendMatrixContour FRET_Matrix_1Bin_mass_conc vs {FRETaxis,size1axis}
SetAxis bottom 0,1
SetAxis left 1,150
ModifyContour FRET_Matrix_1Bin_mass_conc labels=0
ModifyContour FRET_Matrix_1Bin_mass_conc autoLevels={*,*,50}
ModifyContour FRET_Matrix_1Bin_mass_conc ctabLines={0,*,Rainbow256,1}
ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
Label bottom "\Z08 FRET Efficiency"
Label left "\Z08 Size of Oligomer"
//ModifyGraph width=100,height=100
ColorScale/C/N=text4/F=0/A=RC/E frame=0.00
ColorScale/C/N=text4 fsize=8
//ColorScale/C/N=text4 "\Z08Concentration / µM"
TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=0.00/E "\Z08"+y+" hours Mass Distribution (Number of Molecules)"
ModifyGraph manTick(left)={0,20,0,0},manMinor(left)={0,0}
ModifyContour FRET_Matrix_1bin_mass_conc ctabLines={0,0.2,Rainbow256,1},manLevels={0,0.002,100}
ModifyGraph fSize(bottom)=8
ModifyGraph fSize(left)=8
endfor

newLayout/N=Histograms_1
AppendLayoutObject/r=(75,95,280,230)/F=1 graph graph0
AppendLayoutObject/r=(75,229,280,364)/F=1 graph graph2
AppendLayoutObject/r=(75,363,280,498)/F=1 graph graph4
AppendLayoutObject/r=(75,497,280,632)/F=1 graph graph6
AppendLayoutObject/r=(75,631,280,766)/F=1 graph graph8
AppendLayoutObject/r=(279,95,484,230)/F=1 graph graph1
AppendLayoutObject/r=(279,229,484,364)/F=1 graph graph3
AppendLayoutObject/r=(279,363,484,498)/F=1 graph graph5
AppendLayoutObject/r=(279,497,484,632)/F=1 graph graph7
AppendLayoutObject/r=(279,631,484,766)/F=1 graph graph9

newLayout/N=Histograms_2
AppendLayoutObject/r=(75,95,280,230)/F=1 graph graph10
AppendLayoutObject/r=(75,229,280,364)/F=1 graph graph12
AppendLayoutObject/r=(75,363,280,498)/F=1 graph graph14
AppendLayoutObject/r=(75,497,280,632)/F=1 graph graph16
AppendLayoutObject/r=(75,631,280,766)/F=1 graph graph18
AppendLayoutObject/r=(279,95,484,230)/F=1 graph graph11
AppendLayoutObject/r=(279,229,484,364)/F=1 graph graph13
AppendLayoutObject/r=(279,363,484,498)/F=1 graph graph15
AppendLayoutObject/r=(279,497,484,632)/F=1 graph graph17
AppendLayoutObject/r=(279,631,484,766)/F=1 graph graph19


//AppendLayoutObject graph graph2
end


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Deletes last time-point. 

function error()
setdatafolder root:
wave counter
counter[0]=counter[0]-1
end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Counter- housekeeping. 

function count(names)
variable names

SetDataFolder root:
wave timepoints
wave counter


counter[0]+=1
variable n = counter[0]-1
Timepoints[n]=names
end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function sizeplot()
setdatafolder root:
wave timepoints
variable n
for(n=0;n<(dimsize(timepoints,0));n+=1)
setdatafolder root:
variable l=timepoints[n]
string y=num2str(l)
setdatafolder $y

wave data
make/o/n=(dimsize(data,0)) sized=NAN
variable q
for(q=0;q<(dimsize(data,0));q+=1)
sized[q]=data[q][7]
endfor

Make/N=100/O sized_Hist;DelayUpdate
Histogram/P/C/B={0,1,100} sized,sized_Hist
Display/K=0 sized_Hist
SetAxis bottom 0,20
Label left "Probability Density"
SetAxis left 0,0.6
TextBox/C/N=text1/F=0/A=MT/E y+" Hours"
ModifyGraph mode=5,hbFill=3,useNegPat=1,hBarNegFill=3,rgb=(8704,8704,8704)

endfor

newLayout/N=Histograms_1
AppendLayoutObject/r=(75,95,280,230)/F=1 graph graph0
AppendLayoutObject/r=(75,229,280,364)/F=1 graph graph2
AppendLayoutObject/r=(75,363,280,498)/F=1 graph graph4
AppendLayoutObject/r=(75,497,280,632)/F=1 graph graph6
AppendLayoutObject/r=(75,631,280,766)/F=1 graph graph8
AppendLayoutObject/r=(279,95,484,230)/F=1 graph graph1
AppendLayoutObject/r=(279,229,484,364)/F=1 graph graph3
AppendLayoutObject/r=(279,363,484,498)/F=1 graph graph5
AppendLayoutObject/r=(279,497,484,632)/F=1 graph graph7
AppendLayoutObject/r=(279,631,484,766)/F=1 graph graph9

newLayout/N=Histograms_2
AppendLayoutObject/r=(75,95,280,230)/F=1 graph graph10
AppendLayoutObject/r=(75,229,280,364)/F=1 graph graph12
AppendLayoutObject/r=(75,363,280,498)/F=1 graph graph14
AppendLayoutObject/r=(75,497,280,632)/F=1 graph graph16
AppendLayoutObject/r=(75,631,280,766)/F=1 graph graph18
AppendLayoutObject/r=(279,95,484,230)/F=1 graph graph11
AppendLayoutObject/r=(279,229,484,364)/F=1 graph graph13
AppendLayoutObject/r=(279,363,484,498)/F=1 graph graph15
AppendLayoutObject/r=(279,497,484,632)/F=1 graph graph17
AppendLayoutObject/r=(279,631,484,766)/F=1 graph graph19


//AppendLayoutObject graph graph2
end

function allfreter()
setdatafolder root:

wave timepoints
variable n
for(n=0;n<(dimsize(timepoints,0));n+=1)
setdatafolder root:
variable l=timepoints[n]
string y=num2str(l)
setdatafolder $y
duplicate root:AllFret,allFret
wave data,counter1,FRET_Efficiencies_Hist
variable m=counter1[0]
variable p
for(p=0;p<20;p+=1)
allfret[p][m]=FRET_Efficiencies_Hist[p]
endfor
duplicate/o AllFret,root:allFret
endfor

newLayout/N=Histograms_1
AppendLayoutObject/r=(75,95,280,230)/F=1 graph graph0
AppendLayoutObject/r=(75,229,280,364)/F=1 graph graph2
AppendLayoutObject/r=(75,363,280,498)/F=1 graph graph4
AppendLayoutObject/r=(75,497,280,632)/F=1 graph graph6
AppendLayoutObject/r=(75,631,280,766)/F=1 graph graph8
AppendLayoutObject/r=(279,95,484,230)/F=1 graph graph1
AppendLayoutObject/r=(279,229,484,364)/F=1 graph graph3
AppendLayoutObject/r=(279,363,484,498)/F=1 graph graph5
AppendLayoutObject/r=(279,497,484,632)/F=1 graph graph7
AppendLayoutObject/r=(279,631,484,766)/F=1 graph graph9

newLayout/N=Histograms_2
AppendLayoutObject/r=(75,95,280,230)/F=1 graph graph10
AppendLayoutObject/r=(75,229,280,364)/F=1 graph graph12
AppendLayoutObject/r=(75,363,280,498)/F=1 graph graph14
AppendLayoutObject/r=(75,497,280,632)/F=1 graph graph16
AppendLayoutObject/r=(75,631,280,766)/F=1 graph graph18
AppendLayoutObject/r=(279,95,484,230)/F=1 graph graph11
AppendLayoutObject/r=(279,229,484,364)/F=1 graph graph13
AppendLayoutObject/r=(279,363,484,498)/F=1 graph graph15
AppendLayoutObject/r=(279,497,484,632)/F=1 graph graph17
AppendLayoutObject/r=(279,631,484,766)/F=1 graph graph19


//AppendLayoutObject graph graph2
end

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Monomer count

function donorbursts()
setdatafolder root:
make/o/n=(dimsize(timepoints,0)) donorcount
wave timepoints
variable n
for(n=0;n<(dimsize(timepoints,0));n+=1)
variable l = timepoints[n]
string y=num2str(l)
setdatafolder $y
wave fret_donor,donor_bursts
variable monomer=donor_bursts[0]
setdatafolder root:
donorcount[n]=monomer
endfor
end

function change()
wave variables
wave timepoints
variable n
for(n=0;n<(dimsize(timepoints,0));n+=1)
variable l = timepoints[n]
string y=num2str(l)
setdatafolder root:$y
killwaves FRET_Smallmers_Hist
killwaves FRET_Mediummers_Hist
killwaves FRET_Largemers_Hist
killwaves lnz_Smallmers_Hist
killwaves lnz_Mediummers_Hist
killwaves lnz_Largemers_Hist
copy2()
fret_histograms()
lnz_histograms()
endfor
end

function all_FRET()
setdatafolder root:
make/o/n=(20,dimsize(timepoints,0)+1) allFRET

variable m
variable r = 0.025
for(m=0;m<=19;m+=1)
allfret[m][0]=r
r+=0.05
endfor
wave timepoints
variable n
for(n=0;n<(dimsize(timepoints,0));n+=1)
variable l = timepoints[n]
string y=num2str(l)
setdatafolder $y
wave fret_donor
duplicate/o root:allFRET,allFRET
variable j
wave FRET_Efficiencies_hist
for(j=0;j<20;j+=1)
allFRET[j][n+1]=FRET_Efficiencies_hist[j]
endfor
duplicate/o allFRET,root:allFRET

setdatafolder root:

endfor
end

//Consecutive bins test

function consecu()
//variable pod
wave bin_number,file_number,Data
make/o/n=(dimsize(bin_number,0)) consecutive=NAN,tokeep=NAN
make/o/n=(dimsize(bin_number,0),4) keep=NAN,discard=NAN
duplicate/o bin_number,bin_number_edit
make/o/n=(dimsize(data,0),dimsize(data,1)) Data_removed=NAN
make/o/n=(1,(dimsize(data,1))) data_1,data_edit,data_2,data_edit1,data_3,data_4,data_5,data_6,data_7,data_edit2,data_edit3,data_edit4,data_edit5,data_edit6

variable r,s,t
for(s=0;s<(dimsize(data,0));s+=1)
data[s][0]=s
endfor

variable a,b,c,d
for(a=0;a<(dimsize(bin_number,0));a+=1)
consecutive[a]=bin_number[a+1]-bin_number[a]
if(consecutive[a]==1)
//discard[d][0]=file_number[a]
//discard[d][1]=bin_number[a]
//discard[d][2]=a
//discard[d+1][0]=file_number[a+1]
//discard[d+1][1]=bin_number[a+1]
//discard[d+1][2]=a+1
data[a][4]=2
data[a+1][4]=2
d+=2
endif
endfor
b=0
c=0
for(a=0;a<(dimsize(data,0));a+=1)
if(data[a][4]<2)
redimension/n=((b+1),-1) data_1
data_1[b][0]=data[a][0]
data_1[b][1]=data[a][1]
data_1[b][2]=data[a][2]
data_1[b][3]=data[a][3]
data_1[b][4]=data[a][4]
data_1[b][5]=data[a][5]
data_1[b][6]=data[a][6]
data_1[b][7]=data[a][7]
data_1[b][8]=data[a][8]
b+=1
else
redimension/n=((c+1),-1) data_edit
data_edit[c][0]=data[a][0]
data_edit[c][1]=data[a][1]
data_edit[c][2]=data[a][2]
data_edit[c][3]=data[a][3]
data_edit[c][4]=data[a][4]
data_edit[c][5]=data[a][5]
data_edit[c][6]=data[a][6]
data_edit[c][7]=data[a][7]
data_edit[c][8]=data[a][8]
c+=1
endif
endfor
variable e=0
for(a=0;a<(dimsize(data_edit,0));a+=1)
if(data_edit[a+2][5]-data_edit[a][5]==2)
data_edit[a][4]=3
data_edit[a+1][4]=3
data_edit[a+2][4]=3
endif
endfor
b=0
c=0 
for(a=0;a<(dimsize(data_edit,0));a+=1)
if(data_edit[a][4]<3)
redimension/n=((b+1),-1) data_2
data_2[b][0]=data_edit[a][0]
data_2[b][1]=data_edit[a][1]
data_2[b][2]=data_edit[a][2]
data_2[b][3]=data_edit[a][3]
data_2[b][4]=data_edit[a][4]
data_2[b][5]=data_edit[a][5]
data_2[b][6]=data_edit[a][6]
data_2[b][7]=data_edit[a][7]
data_2[b][8]=data_edit[a][8]
b+=1
else
redimension/n=((c+1),-1) data_edit1
data_edit1[c][0]=data_edit[a][0]
data_edit1[c][1]=data_edit[a][1]
data_edit1[c][2]=data_edit[a][2]
data_edit1[c][3]=data_edit[a][3]
data_edit1[c][4]=data_edit[a][4]
data_edit1[c][5]=data_edit[a][5]
data_edit1[c][6]=data_edit[a][6]
data_edit1[c][7]=data_edit[a][7]
data_edit1[c][8]=data_edit[a][8]
c+=1
endif
endfor

make/o/n=(dimsize(data_2,0)/2,9) data_2bin
b=0
for(a=0;a<(dimsize(data_2,0));a+=2)

data_2bin[b][0]=b
data_2bin[b][1]=data_2[a][1]+data_2[a+1][1] //donor
data_2bin[b][2]=data_2[a][2]+data_2[a+1][2] //acceptor
data_2bin[b][4]=data_2[a][4]
b+=1
endfor
wave gamma_factor,mb1
variable mb=mb1[0]
variable gammafact=gamma_factor[0]
for(a=0;a<(dimsize(data_2bin,0));a+=1)
	Data_2bin[a][6] = data_2bin[a][2]/(gammafact*data_2bin[a][1]+data_2bin[a][2]) // FRET Efficiencies
	Data_2bin[a][7] = (2*((data_2bin[a][1]+(data_2bin[a][2]/gammafact))/mb)) 
	Data_2bin[a][8] = ln(data_2bin[a][2]/data_2bin[a][1])
endfor

e=0
for(a=0;a<(dimsize(data_edit1,0));a+=1)
if(data_edit1[a+3][5]-data_edit1[a][5]==3)
data_edit1[a][4]=4
data_edit1[a+1][4]=4
data_edit1[a+2][4]=4
data_edit1[a+3][4]=4
endif
endfor
b=0
c=0

for(a=0;a<(dimsize(data_edit1,0));a+=1)
if(data_edit1[a][4]<4)
redimension/n=((b+1),-1) data_3
data_3[b][0]=data_edit1[a][0]
data_3[b][1]=data_edit1[a][1]
data_3[b][2]=data_edit1[a][2]
data_3[b][3]=data_edit1[a][3]
data_3[b][4]=data_edit1[a][4]
data_3[b][5]=data_edit1[a][5]
data_3[b][6]=data_edit1[a][6]
data_3[b][7]=data_edit1[a][7]
data_3[b][8]=data_edit1[a][8]
b+=1
else
redimension/n=((c+1),-1) data_edit2
data_edit2[c][0]=data_edit1[a][0]
data_edit2[c][1]=data_edit1[a][1]
data_edit2[c][2]=data_edit1[a][2]
data_edit2[c][3]=data_edit1[a][3]
data_edit2[c][4]=data_edit1[a][4]
data_edit2[c][5]=data_edit1[a][5]
data_edit2[c][6]=data_edit1[a][6]
data_edit2[c][7]=data_edit1[a][7]
data_edit2[c][8]=data_edit1[a][8]
c+=1
endif
endfor

make/o/n=(dimsize(data_3,0)/3,9) data_3bin
b=0
for(a=0;a<(dimsize(data_3,0));a+=3)

data_3bin[b][0]=b
data_3bin[b][1]=data_3[a][1]+data_3[a+1][1]+data_3[a+2][1] //donor
data_3bin[b][2]=data_3[a][2]+data_3[a+1][2]+data_3[a+2][2] //acceptor
data_3bin[b][4]=data_3[a][4]
b+=1
endfor
for(a=0;a<(dimsize(data_3bin,0));a+=1)
	Data_3bin[a][6] = data_3bin[a][2]/(gammafact*data_3bin[a][1]+data_3bin[a][2]) // FRET Efficiencies
	Data_3bin[a][7] = (2*((data_3bin[a][1]+(data_3bin[a][2]/gammafact))/mb)) 
	Data_3bin[a][8] = ln(data_3bin[a][2]/data_3bin[a][1])
endfor

e=0
for(a=0;a<(dimsize(data_edit2,0));a+=1)
if(data_edit2[a+4][5]-data_edit2[a][5]==4)
data_edit2[a][4]=5
data_edit2[a+1][4]=5
data_edit2[a+2][4]=5
data_edit2[a+3][4]=5
data_edit2[a+4][4]=5
endif
endfor

b=0
c=0

for(a=0;a<(dimsize(data_edit2,0));a+=1)
if(data_edit2[a][4]<5)
redimension/n=((b+1),-1) data_4
data_4[b][0]=data_edit2[a][0]
data_4[b][1]=data_edit2[a][1]
data_4[b][2]=data_edit2[a][2]
data_4[b][3]=data_edit2[a][3]
data_4[b][4]=data_edit2[a][4]
data_4[b][5]=data_edit2[a][5]
data_4[b][6]=data_edit2[a][6]
data_4[b][7]=data_edit2[a][7]
data_4[b][8]=data_edit2[a][8]
b+=1
else
redimension/n=((c+1),-1) data_edit3
data_edit3[c][0]=data_edit2[a][0]
data_edit3[c][1]=data_edit2[a][1]
data_edit3[c][2]=data_edit2[a][2]
data_edit3[c][3]=data_edit2[a][3]
data_edit3[c][4]=data_edit2[a][4]
data_edit3[c][5]=data_edit2[a][5]
data_edit3[c][6]=data_edit2[a][6]
data_edit3[c][7]=data_edit2[a][7]
data_edit3[c][8]=data_edit2[a][8]
c+=1
endif
endfor

b=0
make/o/n=(dimsize(data_4,0)/4,9) data_4bin
for(a=0;a<(dimsize(data_4,0));a+=4)

data_4bin[b][0]=b
data_4bin[b][1]=data_4[a][1]+data_4[a+1][1]+data_4[a+2][1]+data_4[a+3][1] //donor
data_4bin[b][2]=data_4[a][2]+data_4[a+1][2]+data_4[a+2][2]+data_4[a+3][2] //acceptor
data_4bin[b][4]=data_4[a][4]
b+=1
endfor

for(a=0;a<(dimsize(data_4bin,0));a+=1)
	Data_4bin[a][6] = data_4bin[a][2]/(gammafact*data_4bin[a][1]+data_4bin[a][2]) // FRET Efficiencies
	Data_4bin[a][7] = (2*((data_4bin[a][1]+(data_4bin[a][2]/gammafact))/mb)) 
	Data_4bin[a][8] = ln(data_4bin[a][2]/data_4bin[a][1])
endfor

e=0
for(a=0;a<(dimsize(data_edit3,0));a+=1)
if(data_edit3[a+4][5]-data_edit3[a][5]==5)
data_edit3[a][4]=6
data_edit3[a+1][4]=6
data_edit3[a+2][4]=6
data_edit3[a+3][4]=6
data_edit3[a+4][4]=6
endif
endfor

b=0
c=0

for(a=0;a<(dimsize(data_edit3,0));a+=1)
if(data_edit3[a][4]<6)
redimension/n=((b+1),-1) data_5
data_5[b][0]=data_edit3[a][0]
data_5[b][1]=data_edit3[a][1]
data_5[b][2]=data_edit3[a][2]
data_5[b][3]=data_edit3[a][3]
data_5[b][4]=data_edit3[a][4]
data_5[b][5]=data_edit3[a][5]
data_5[b][6]=data_edit3[a][6]
data_5[b][7]=data_edit3[a][7]
data_5[b][8]=data_edit3[a][8]
b+=1
else
redimension/n=((c+1),-1) data_edit4
data_edit4[c][0]=data_edit3[a][0]
data_edit4[c][1]=data_edit3[a][1]
data_edit4[c][2]=data_edit3[a][2]
data_edit4[c][3]=data_edit3[a][3]
data_edit4[c][4]=data_edit3[a][4]
data_edit4[c][5]=data_edit3[a][5]
data_edit4[c][6]=data_edit3[a][6]
data_edit4[c][7]=data_edit3[a][7]
data_edit4[c][8]=data_edit3[a][8]
c+=1
endif
endfor

b=0
make/o/n=(dimsize(data_5,0)/5,9) data_5bin
for(a=0;a<(dimsize(data_5,0));a+=5)

data_5bin[b][0]=b
data_5bin[b][1]=data_5[a][1]+data_5[a+1][1]+data_5[a+2][1]+data_5[a+3][1]+data_5[a+4][1] //donor
data_5bin[b][2]=data_5[a][2]+data_5[a+1][2]+data_5[a+2][2]+data_5[a+3][2]+data_5[a+4][2] //acceptor
data_5bin[b][4]=data_5[a][4]
b+=1
endfor

for(a=0;a<(dimsize(data_5bin,0));a+=1)
	Data_5bin[a][6] = data_5bin[a][2]/(gammafact*data_5bin[a][1]+data_5bin[a][2]) // FRET Efficiencies
	Data_5bin[a][7] = (2*((data_5bin[a][1]+(data_5bin[a][2]/gammafact))/mb)) 
	Data_5bin[a][8] = ln(data_5bin[a][2]/data_5bin[a][1])
endfor

make/o/n=(1,9) data_all

variable full=0

for(a=0;a<(dimsize(data_1,0));a+=1)
redimension/n=(full+1,-1) data_all
data_all[full][0]=data_1[a][0]
data_all[full][1]=data_1[a][1]
data_all[full][2]=data_1[a][2]
data_all[full][3]=data_1[a][3]
data_all[full][4]=data_1[a][4]
data_all[full][5]=data_1[a][5]
data_all[full][6]=data_1[a][6]
data_all[full][7]=data_1[a][7]
data_all[full][8]=data_1[a][8]
full+=1
endfor
// The following code adds the wider events to the dataset which is to be analysed. If this is not needed, then
// comment out.

//for(a=0;a<(dimsize(data_2bin,0));a+=1)
//redimension/n=(full+1,-1) data_all
//data_all[full][0]=data_2bin[a][0]
//data_all[full][1]=data_2bin[a][1]
//data_all[full][2]=data_2bin[a][2]
//data_all[full][3]=data_2bin[a][3]
//data_all[full][4]=data_2bin[a][4]
//data_all[full][5]=data_2bin[a][5]
//data_all[full][6]=data_2bin[a][6]
//data_all[full][7]=data_2bin[a][7]
//data_all[full][8]=data_2bin[a][8]
//full+=1
//endfor

//for(a=0;a<(dimsize(data_3bin,0));a+=1)
//redimension/n=(full+1,-1) data_all
//data_all[full][0]=data_3bin[a][0]
//data_all[full][1]=data_3bin[a][1]
//data_all[full][2]=data_3bin[a][2]
//data_all[full][3]=data_3bin[a][3]
///data_all[full][4]=data_3bin[a][4]
//data_all[full][5]=data_3bin[a][5]
//data_all[full][6]=data_3bin[a][6]
//data_all[full][7]=data_3bin[a][7]
//data_all[full][8]=data_3bin[a][8]
//full+=1
//endfor

//for(a=0;a<(dimsize(data_4bin,0));a+=1)
//redimension/n=(full+1,-1) data_all
//data_all[full][0]=data_4bin[a][0]
//data_all[full][1]=data_4bin[a][1]
//data_all[full][2]=data_4bin[a][2]
//data_all[full][3]=data_4bin[a][3]
//data_all[full][4]=data_4bin[a][4]
//data_all[full][5]=data_4bin[a][5]
//data_all[full][6]=data_4bin[a][6]
//data_all[full][7]=data_4bin[a][7]
//data_all[full][8]=data_4bin[a][8]
//full+=1
//endfor

// End of adding wider bins.

duplicate/o data_all,data
end




function normalise()
wave timepoints,specA,specB,small,mediumhigh,mediumlow,large,all
duplicate/o specA,specAN
duplicate/o specB,specBN
duplicate/o small,smallN
duplicate/o mediumlow,mediumlowN
duplicate/o mediumhigh,mediumhighN
duplicate/o large,largeN
duplicate/o all,alln
donorbursts()
wave donorcount
variable a
for(a=0;a<(dimsize(timepoints,0));a+=1)
specAN[a]=specan[a]/donorcount[a]
specbN[a]=specbn[a]/donorcount[a]
smallN[a]=smalln[a]/donorcount[a]
mediumlowN[a]=mediumlown[a]/donorcount[a]
mediumhighN[a]=mediumhighn[a]/donorcount[a]
largeN[a]=largen[a]/donorcount[a]
allN[a]=alln[a]/donorcount[a]
endfor
Display smallN,mediumlowN,mediumhighN,largeN,alln vs timepoints
ModifyGraph notation(left)=1;DelayUpdate
Label left "Fraction of Oligomers";DelayUpdate
Label bottom "Time / hours"
ModifyGraph mode=4,marker=19,msize=1.5,rgb(smallN)=(65280,49152,16384);DelayUpdate
ModifyGraph rgb(mediumlowN)=(0,52224,0),rgb(mediumhighN)=(39168,0,31232);DelayUpdate
ModifyGraph rgb(largeN)=(39168,0,0)
Legend/C/N=text0/J/F=0/A=RC/E "\\s(smallN) Smallmers\r\\s(mediumlowN) Mediummers (low FRET)\r\\s(mediumhighN) Mediummers (high FRET)\r\\s(largeN) Largemers";DelayUpdate
AppendText "\\s(alln) All"


Display specAN,specBN vs timepoints
Label bottom "Time / hours"
Label left "Fraction of oligomers";DelayUpdate
SetAxis left 0,*
ModifyGraph notation(left)=1
ModifyGraph mode=4,marker=19,msize=1.5,rgb(specAN)=(0,0,65280)
Legend/C/N=text0/J/F=0/A=RC/E "\\s(specAN) Species A\r\\s(specBN) Species B"
end

Function add(hours,file,donorthresh,acceptorthresh,first,last,bins,binwidth)
variable hours,donorthresh,acceptorthresh,first,last,bins,binwidth
string file
variable m,n
setdatafolder root:
wave timepoints,variables,lnz_large,lnz_medium,lnz_small,verylarge,FRET_small,FRET_medium,FRET_large,counter,width_deletion,allfret
redimension/n=(-1,dimsize(fret_small,1)+1)  allfret,lnz_large,lnz_medium,lnz_small,FRET_small,FRET_medium,FRET_large,verylarge
variable add_numbers=(dimsize(fret_small,1)+1)
redimension/n=(dimsize(timepoints,0)+1) timepoints,width_deletion
timepoints[add_numbers]=hours
string nameit=num2str(hours)
wave counter
counter[0]+=1
newdatafolder/s/o  $nameit
copy()
make/o/n=1 hour=hours

wave variables
variable pod=(dimsize(timepoints,0)+1)
loader2(file,donorthresh,acceptorthresh,first,last,bins,binwidth)

get()
consecu()
zplot()
matrices_fret()
FRET_histograms()
//matrices_z()
//znz_histograms()
SetDataFolder root:

end






////////////////////////////7


function loader2(filename,donorthresh,acceptorthresh,first,last,bins,binwidth)

variable donorthresh
variable acceptorthresh
variable first
variable last
variable bins
variable binwidth
string filename
variable autodonor=0.672
variable autoacceptor=0.541
variable crosstalk=0.080
variable gammafactor=0.26



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
variable m,n=0
variable d=bins*(last-first)
variable frametime=(binwidth/1000000000)*bins
make/o/n=1 gamma_factor=gammafactor
print frametime
variable g,h
make/o/n=(1) FRET_donor=NAN,FRET_acceptor=NAN,donorchannel=NAN,bin_number=NAN
make/o/n=((last-first),3) Stats
LoadWave/O/J/D/W/K=0 	// Loads wave to get path
string path=S_path+filename
print path
string num
variable c
for(c=first;c<=last;c+=1)
if(c<10)
num="000"+num2str(c)
elseif(c<100)
num="00"+num2str(c)
elseif(c<1000)
num="0"+num2str(c)
endif

string path2=path+num+".dat"

GBLoadWave/V/T={32,4}/W=2 path2
wave wave0,wave1,wave2
variable e,k,f
for(e=0;e<=(dimsize(wave0,0));e+=1)
if(wave0[e]>donorthresh&&acceptorthresh<(wave1[e]-crosstalk*wave0[e]))
Redimension/N=(k) FRET_Donor,FRET_Acceptor,bin_number
FRET_Donor[k]=wave0[e]-autodonor
FRET_Acceptor[k]=wave1[e]-autoacceptor-crosstalk*wave0[e]
bin_number[k]=e
k+=1
m+=1
elseif(wave0[e]>donorthresh)
Redimension/N=(f) donorchannel
donorchannel[f]=wave0[e]
f+=1
n+=1
endif


endfor
stats[c][0]=n/frametime
stats[c][1]=m/frametime
stats[c][2]=m/n
m=0
n=0

killwaves wave0
killwaves wave1
killwaves wave2

endfor

display stats[][0]
TextBox/C/N=text0/F=0/A=MT/E Pathway
Label bottom "Frame Number"
Label left "Event rate (s\\S-1\\M)"
AppendToGraph/R stats[][1]
Label right "Event rate (s\\S-1\\M)"
ModifyGraph rgb(Stats#1)=(0,0,65280)
ModifyGraph axRGB(right)=(0,0,65280)
ModifyGraph axRGB(left)=(65280,3328,0),tlblRGB(left)=(65280,3328,0);DelayUpdate
ModifyGraph tlblRGB(right)=(0,0,65280),alblRGB(left)=(65280,3328,0);DelayUpdate
ModifyGraph alblRGB(right)=(0,0,65280)
make/o/n=1 donor_bursts=k+f
wavestats donorchannel
make/o/n=1 monomerb=V_avg
SetAxis left 0,*
SetAxis right 0,*
Legend/C/N=text1/J/F=2/A=RC/E "\\s(Stats) Donor\r\\s(Stats#1) Coincident"
print V_avg

end


function allfretfit()
setdatafolder root:
wave allfret,timepoints
variable number=(dimsize(allfret,1))
variable a
make/o/n=(7,(dimsize(timepoints,0))) coeffs
make/o/n=(dimsize(timepoints,0)) allLow,allHigh,allLow_norm,allHigh_norm
for(a=1;a<(number);a+=1)
display allfret[][a] vs allfret[][0]
ModifyGraph mode=5,hbFill=4,useNegPat=1,hBarNegFill=4,rgb=(0,0,0)
Label bottom "FRET Efficiency"
Label left "Number of Events"

string fit="fit"+num2str(a)
make/o/n=20 $fit
coeffs[][a-1]={0,30,0.235,0.065,50,0.4196,0.09021}
FuncFit/H="1011011"/NTHR=0/TBOX=768 DoubleGaussian coeffs[][a-1]  AllFret[*][a]/X=AllFret[*][0]/D=$fit
//FuncFit/NTHR=0/TBOX=768 DoubleGaussian coeffs[][a-1]  AllFret[*][a]/X=AllFret[*][0]/D=$fit
appendtograph $fit vs allfret[][0]
allLow[a-1]=coeffs[1][a-1]*(sqrt(pi/coeffs[3][a-1]))
allhigh[a-1]=coeffs[4][a-1]*(sqrt(pi/coeffs[6][a-1]))
wave donorcount
allLow_norm[a-1]=coeffs[1][a-1]*(sqrt(pi/coeffs[3][a-1]))/donorcount[a-1]
allhigh_norm[a-1]=coeffs[4][a-1]*(sqrt(pi/coeffs[6][a-1]))/donorcount[a-1]


endfor
Display allHigh,allLow vs timepoints
ModifyGraph mode=4,marker=19,msize=1.5,rgb(allHigh)=(65280,0,0);DelayUpdate
ModifyGraph rgb(allLow)=(0,15872,65280)
Label bottom "Time / hours"
Label left "Number of events"

end


function size(FRET_cutoff)
variable fret_cutoff
setdatafolder root:
variable d
string name
wave timepoints
make/o/n=(50,dimsize(timepoints,0)+1) sizes_high=0,sizes_low=0

variable n
for(n=0;n<(dimsize(timepoints,0));n+=1)
variable t=1
t+=1

variable l = timepoints[n]
string y=num2str(l)
setdatafolder root:
setdatafolder $y
copyit()
wave data
variable x,z
variable f=1
	for(d=0;d<50;d+=1)
		variable count=0
		for(z=0;z<(dimsize(data,0));z+=1)
			if(data[z][7]<(d+2) && data[z][6]<FRET_cutoff)
			if(data[z][7]>d)
			count+=1
			endif
			endif
			sizes_low[d][0]=d
			
		endfor
		wave sizes_low
sizes_low[f][n]=count

count=0
f+=1
endfor
count=0
f=1
	for(d=0;d<50;d+=1)
		for(z=0;z<(dimsize(data,0));z+=1)
			if(data[z][7]<(d+2) && data[z][6]>FRET_cutoff)
			if(data[z][7]>d && data[z][6]>FRET_cutoff)
			count+=1
			endif
			endif
		sizes_high[d][0]=d
		endfor
		
		wave sizes_high
sizes_high[f][n]=count
//print count
count=0
f+=1
endfor

copyit2()

endfor

setdatafolder root:
end

function copyit()
wave sizes_low,sizes_high
duplicate/o root:sizes_low,sizes_low
duplicate/o root:sizes_high,sizes_high
end
function copyit2()
wave sizes_low,sizes_high
duplicate/o sizes_low,root:sizes_low
duplicate/o sizes_high,root:sizes_high

end




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Function to see size distribution /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function seesize(times)
variable times		
kill()				// Kill any open windows
setdatafolder root:	// Set data folder to root
wave sizes_low,sizes_high		// Let program know which waves will be looked at. 

display sizes_low[][times] vs sizes_low[][0]	// show nice graph of sizes. 
// Make pretty. 
appendtograph sizes_high[][times] vs sizes_high[][0]
Label bottom "Approximate Size"
Label left "Number of counts"
ModifyGraph lsize=2,rgb(sizes_low)=(0,0,26112)
Legend/C/N=text0/J/F=0/A=RC/E "\\s(sizes_low) Low FRET\r\\s(sizes_high) High FRET"
ModifyGraph width=212.598,height={Aspect,0.75},gfSize=15
SetAxis left 0,350

end


Function seeit()
	setdatafolder root:
	wave timepoints
	String traceName
	string Hours
	variable length=dimsize(timepoints,0)
	if(length<1)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])
		
	elseif(length<2)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])
			
		elseif(length<3)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])
		
		elseif(length<4)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])+";"+num2str(timepoints[2])
			
			elseif(length<5)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])+";"+num2str(timepoints[2])+";"+num2str(timepoints[3])
		
			elseif(length<6)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])+";"+num2str(timepoints[2])+";"+num2str(timepoints[3])+";"+num2str(timepoints[4])
			elseif(length<7)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])+";"+num2str(timepoints[2])+";"+num2str(timepoints[3])+";"+num2str(timepoints[4])+";"+num2str(timepoints[5])
			elseif(length<8)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])+";"+num2str(timepoints[2])+";"+num2str(timepoints[3])+";"+num2str(timepoints[4])+";"+num2str(timepoints[5])+";"+num2str(timepoints[6])
		elseif(length<9)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])+";"+num2str(timepoints[2])+";"+num2str(timepoints[3])+";"+num2str(timepoints[4])+";"+num2str(timepoints[5])+";"+num2str(timepoints[6])+";"+num2str(timepoints[7])
		elseif(length<10)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])+";"+num2str(timepoints[2])+";"+num2str(timepoints[3])+";"+num2str(timepoints[4])+";"+num2str(timepoints[5])+";"+num2str(timepoints[6])+";"+num2str(timepoints[7])+";"+num2str(timepoints[8])
		elseif(length<11)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])+";"+num2str(timepoints[2])+";"+num2str(timepoints[3])+";"+num2str(timepoints[4])+";"+num2str(timepoints[5])+";"+num2str(timepoints[6])+";"+num2str(timepoints[7])+";"+num2str(timepoints[8])+";"+num2str(timepoints[9])
		elseif(length<12)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])+";"+num2str(timepoints[2])+";"+num2str(timepoints[3])+";"+num2str(timepoints[4])+";"+num2str(timepoints[5])+";"+num2str(timepoints[6])+";"+num2str(timepoints[7])+";"+num2str(timepoints[8])+";"+num2str(timepoints[9])+";"+num2str(timepoints[10])
		elseif(length<13)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])+";"+num2str(timepoints[2])+";"+num2str(timepoints[3])+";"+num2str(timepoints[4])+";"+num2str(timepoints[5])+";"+num2str(timepoints[6])+";"+num2str(timepoints[7])+";"+num2str(timepoints[8])+";"+num2str(timepoints[9])+";"+num2str(timepoints[10])+";"+num2str(timepoints[11])
		elseif(length<14)
		Prompt hours,"Timepoint to see:",popup,num2str(timepoints[0])+";"+num2str(timepoints[1])+";"+num2str(timepoints[2])+";"+num2str(timepoints[3])+";"+num2str(timepoints[4])+";"+num2str(timepoints[5])+";"+num2str(timepoints[6])+";"+num2str(timepoints[7])+";"+num2str(timepoints[8])+";"+num2str(timepoints[9])+";"+num2str(timepoints[10])+";"+num2str(timepoints[11])+";"+num2str(timepoints[12])
		elseif(length<15)
		
		
	
	else
	Prompt hours,"Timepoint to see:"
	endif
	DoPrompt "Timepoint",Hours
	
	print hours
	setdatafolder root:$hours
kill()
see()
	
End


/////////////////////////////////////////////   MACROS GO HERE /////////////////////////////////////////////////////////
// Main macro for loading experiments. 
Macro Run(filename,donor,acceptor,first,last,bins,binwidth,x,y)
string filename="Asyn"
Prompt Filename, "Filename (without suffix): "			
variable donor=10
Prompt donor, "Donor Threshold: "			
variable Acceptor=10
Prompt Acceptor, "Acceptor Threshold: "			
variable first=0
Prompt first, "First Filenumber: "			
variable last=79
Prompt last, "Last Filenumber: "			
variable bins=100000
Prompt bins, "Bins in one frame: "			
variable binwidth=50000
Prompt binwidth, "Bin Width (ns): "			
variable x=3
Prompt x, "Set upper limit of first size distribution: "			
Variable y=150
Prompt y, "Set upper limit of second size distribution: "		


experiment(filename,donor,acceptor,first,last,bins,x,y,binwidth)
end

// Macro to show graphs

macro Kinetics()
totalsfret()
normalise()
end

// Macro to show 2D plots

macro Plots()
twod()
end

// Macro shows sizes at different FRET efficiencies.

Macro Size_Split(FRET)
Variable FRET=0.6
Prompt FRET, "FRET cutoff: "			// Set prompt for X param


size(FRET)
end

// Macro shows nice layout for plots

macro Timepoint()
seeit()
end

// Macro shows size distribution at time-points

macro See_size(times)
Variable times=0
Prompt times, "Timepoint: "			// Set prompt for X param


seesize(times)
end

// Macro to change boundaries of the FRET histograms.

macro Change_Limits(lower,upper)
variable lower
prompt lower,"New upper limit of first size distribution"
variable upper
prompt upper,"New upper limit of second size distribution"
setdatafolder root:
change_limitsf(lower,upper)
change()
end

function change_limitsf(lower,upper)
variable upper,lower
setdatafolder root:
wave variables 
variables[1]=lower
variables[2]=upper
end 


Macro Add_Timepoint(hours,filename,donor,acceptor,first,last,bins,binwidth)
Variable hours=0
Prompt hours, "Hours: "		// Set prompt for y param
string filename="Asyn"
Prompt Filename, "Filename (without suffix): "			// Set prompt for X param
variable donor=10
Prompt donor, "Donor Threshold: "			// Set prompt for X param
variable Acceptor=10
Prompt Acceptor, "Acceptor Threshold: "			// Set prompt for X param
variable first=0
Prompt first, "First Filenumber: "			// Set prompt for X param
variable last=100
Prompt last, "Last Filenumber: "			// Set prompt for X param
variable bins=100000
Prompt bins, "Bins in one frame: "			// Set prompt for X param
variable binwidth=100000
Prompt binwidth, "Bin Width (ns): "			// Set prompt for X param

add(hours,filename,donor,acceptor,first,last,bins,binwidth)
end

macro Close_Windows()
kill()
end


function lnzfigure2()
wave data
variable length=(dimsize(data,0))
make/o/n=1 zplotdat
variable a
variable c=1
for(a=0;a<(length);a+=1)
if(data[a][7]>5)
redimension/n=(c) zplotdat
zplotdat[a]=data[a][8]
c+=1
endif
endfor
Make/N=30/O zplotdat_Hist;DelayUpdate
Histogram/B={-3,0.2,30} zplotdat,zplotdat_Hist
make/o/n=30 zploty,zPlotyFRET
variable b=-3
for(a=0;a<30;a+=1)
zploty[a]=b
zplotyFRET[a]=exp(b)/(1+exp(b))
b+=0.2
endfor
Make/D/N=7/O W_coef
W_coef[0] = {0,50,-1,1,150,0,1}
FuncFit/NTHR=0/TBOX=768 Double_Gauss_Amp W_coef  zplotdat_Hist /X=zploty /D 
make/o/n=1 twochi=V_chisq
duplicate/o  fit_zplotdat_Hist,doublefit
make/o/n=200 fit_axis,fit_axis_fret
variable val=-3
for(a=0;a<200;a+=1)
fit_axis[a]=val
fit_axis_fret[a]=exp(val)/(1+exp(val))
val+=0.029

endfor

make/o/n=(4) peak1={w_coef[0],w_coef[1],w_coef[2],w_coef[3]},peak2={w_coef[0],w_coef[4],w_coef[5],w_coef[6]}
FuncFit/H="1111"/NTHR=0/TBOX=768 Gauss_Amp peak1  zplotdat_Hist /X=zploty /D 
duplicate/o  fit_zplotdat_Hist,fit1
FuncFit/H="1111"/NTHR=0/TBOX=768 Gauss_Amp peak2  zplotdat_Hist /X=zploty /D 
duplicate/o  fit_zplotdat_Hist,fit2

display zplotdat_Hist vs zPlotyFRET
AppendToGraph fit1 vs fit_axis_fret
AppendToGraph fit2 vs fit_axis_fret
AppendToGraph doublefit vs fit_axis_fret
ModifyGraph marker=19,rgb=(0,0,39168),mode(zplotdat_Hist)=3;DelayUpdate
Label left "no. of counts";DelayUpdate
Label bottom "FRET Efficiency"
ModifyGraph width=283.465,height=198.425,gFont="Arial",gfSize=15
ModifyGraph rgb(zplotdat_Hist)=(0,0,26112),lsize(fit1)=2,rgb(fit1)=(0,39168,0);DelayUpdate
ModifyGraph lsize(fit2)=2,rgb(fit2)=(65280,21760,0),lsize(doublefit)=2;DelayUpdate
ModifyGraph rgb(doublefit)=(0,0,26112)
wave w_coef
string cent1=num2str(exp(w_coef[2])/(1+exp(w_coef[2])))
string width1=num2str(exp(w_coef[3])/(1+exp(w_coef[3])))
string cent2=num2str(exp(w_coef[5])/(1+exp(w_coef[5])))
string width2=num2str(exp(w_coef[6])/(1+exp(w_coef[6])))
TextBox/C/N=text3/F=0/A=MT/E "Centre = "+cent1+","+cent2+"\rWidth = "+width1+","+width2
end


function lnzfigure1()
wave data
variable length=(dimsize(data,0))
make/o/n=(length) zplotdat
variable c=1
variable a
for(a=0;a<(length);a+=1)
if(data[a][7]>5)
redimension/n=(c) zplotdat
zplotdat[a]=data[a][8]
c+=1
endif
endfor
Make/N=30/O zplotdat_Hist;DelayUpdate
Histogram/B={-3,0.2,30} zplotdat,zplotdat_Hist
make/o/n=30 zploty,zPlotyFRET
variable b=-3
for(a=0;a<30;a+=1)
zploty[a]=b
zplotyFRET[a]=exp(b)/(1+exp(b))
b+=0.2
endfor
CurveFit/NTHR=0/TBOX=768 gauss zplotdat_Hist /D 
make/o/n=1 onechi=V_chisq
make/o/n=200 fit_axis,fit_axis_fret
variable val=-3
for(a=0;a<200;a+=1)
fit_axis[a]=val
fit_axis_fret[a]=exp(val)/(1+exp(val))
val+=0.029
endfor

Display zplotdat_Hist vs zPlotyFRET
AppendToGraph fit_zplotdat_Hist vs fit_zplotdat_Hist
RemoveFromGraph fit_zplotdat_Hist
AppendToGraph fit_zplotdat_Hist vs fit_axis_fret
Label left "no. of counts";DelayUpdate
Label bottom "FRET Efficiency"
ModifyGraph width=283.465,height=198.425,gFont="Arial",gfSize=15
ModifyGraph rgb(zplotdat_Hist)=(0,0,26112),lsize(fit_zplotdat_Hist)=2;DelayUpdate
ModifyGraph rgb(fit_zplotdat_Hist)=(0,39168,0)
ModifyGraph rgb=(0,0,26112),mode(zplotdat_Hist)=3,marker(zplotdat_Hist)=1
wave w_coef
string cent=num2str(exp(w_coef[2])/(1+exp(w_coef[2])))
string width=num2str(exp(w_coef[3])/(1+exp(w_coef[3])))
TextBox/C/N=text3/F=0/A=MT/E "Centre = "+cent+"\rWidth = "+width
end

function lnzfit1()
wave data
variable length=(dimsize(data,0))
make/o/n=(length) zplotdat
variable c=1
variable a
for(a=0;a<(length);a+=1)
if(data[a][7]>5)
redimension/n=(c) zplotdat
zplotdat[a]=data[a][8]
c+=1
endif
endfor
Make/N=30/O zplotdat_Hist;DelayUpdate
Histogram/B={-3,0.2,30} zplotdat,zplotdat_Hist
make/o/n=30 zploty,zPlotyFRET
variable b=-3
for(a=0;a<30;a+=1)
zploty[a]=b
zplotyFRET[a]=exp(b)/(1+exp(b))
b+=0.2
endfor
CurveFit/NTHR=0/TBOX=768 gauss zplotdat_Hist /D 
wavestats zplotdat_hist
make/o/n=1 onechi=V_chisq

end

function lnzfit2()
wave data
variable length=(dimsize(data,0))
make/o/n=1 zplotdat
variable a
variable c=1
for(a=0;a<(length);a+=1)
if(data[a][7]>5)
redimension/n=(c) zplotdat
zplotdat[a]=data[a][8]
c+=1
endif
endfor
Make/N=30/O zplotdat_Hist;DelayUpdate
Histogram/B={-3,0.2,30} zplotdat,zplotdat_Hist
make/o/n=30 zploty,zPlotyFRET
variable b=-3
for(a=0;a<30;a+=1)
zploty[a]=b
zplotyFRET[a]=exp(b)/(1+exp(b))
b+=0.2
endfor
make/o/n=7 W_coef
W_coef[0] = {0,50,-1,1,150,0,1}
FuncFit/NTHR=0/TBOX=768 Double_Gauss_Amp W_coef  zplotdat_Hist /D 
make/o/n=1 twochi=V_chisq

end

function compare()
lnzfit1()
lnzfit2()
wave onechi,twochi
variable percentage_decrease=100*(twochi[0]-onechi[0])/(onechi[0])
print percentage_decrease
make/o/n=1 percent_dec=percentage_decrease
//if(percentage_decrease<-20)
//lnzfigure2()
//else
lnzfigure1()
//endif
string percentage=num2str(percent_dec[0])

make/o/n=1 dec=percentage_decrease
TextBox/C/N=text0/A=MB/E "Change = "+percentage+"%"
TextBox/C/N=text0/F=0
ModifyGraph mode(zplotdat_Hist)=5
ModifyGraph rgb(zplotdat_Hist)=(26112,26112,26112)
SetAxis left 0,*
end

function compare_all()
setdatafolder root:
wave timepoints
variable a
for(a=0;a<(dimsize(timepoints,0));a+=1)
string folder=num2str(timepoints[a])
setdatafolder $folder
compare()
TextBox/C/N=text1 ""+folder+" hours"
setdatafolder root:
endfor
end

Function Gauss_Amp(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0+A*exp(-0.5*((x-xc)/w)^2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = y0
	//CurveFitDialog/ w[1] = A
	//CurveFitDialog/ w[2] = xc
	//CurveFitDialog/ w[3] = w

	return w[0]+w[1]*exp(-0.5*((x-w[2])/w[3])^2)
End

Function Double_Gauss_Amp(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0+A*exp(-0.5*((x-xc)/w)^2)+A2*exp(-0.5*((x-xc2)/w2)^2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 7
	//CurveFitDialog/ w[0] = y0
	//CurveFitDialog/ w[1] = a
	//CurveFitDialog/ w[2] = xc
	//CurveFitDialog/ w[3] = w
	//CurveFitDialog/ w[4] = a2
	//CurveFitDialog/ w[5] = xc2
	//CurveFitDialog/ w[6] = w2

	return w[0]+w[1]*exp(-0.5*((x-w[2])/w[3])^2)+w[4]*exp(-0.5*((x-w[5])/w[6])^2)
End



function extract_param()					//Extract width data after compare_all
setdatafolder root:
variable d
string name
wave timepoints
make/o/n=(50,dimsize(timepoints,0)+1) sizes_high=0,sizes_low=0
make/o/n=(dimsize(timepoints,0)) centre,width
variable n
for(n=0;n<(dimsize(timepoints,0));n+=1)
variable t=1
t+=1

variable l = timepoints[n]
string y=num2str(l)
setdatafolder root:
setdatafolder $y
duplicate/o root:centre,centre2
duplicate/o root:width,width2
wave W_coef
centre2[n]=exp(W_coef[2])/(1+exp(W_coef[2]))
width2[n]=exp(W_coef[3])/(1+exp(W_coef[3]))
duplicate/o centre2,root:centre
duplicate/o width2,root:width




endfor

setdatafolder root:
end