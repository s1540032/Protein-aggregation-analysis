#pragma rtGlobals=1		// Use modern global access method.

Macro Run(foldername,filename,donor,acceptor,last,autod,autob,crosst,crosst2)
string foldername=""
Prompt foldername, "Experiment reference: "
string filename=""
Prompt Filename, "Filename (without suffix): "			// Set prompt for X param
variable donor=5
Prompt donor, "Channel A Threshold: "			// Set prompt for X param
variable Acceptor=5
Prompt Acceptor, "Channel B Threshold: "			// Set prompt for X param
variable last=10
Prompt last, "Number of Files: "			// Set prompt for X param
variable Autod=0.272
Prompt Autod, "Autofluorescence in channel A: "			// Set prompt for X param
variable Autob=1.93
Prompt Autob, "Autofluorescence in channel B: "			// Set prompt for X param
variable crosst=0.099
Prompt crosst, "Blue --> Red crosstalk"			// Set prompt for X param
variable crosst2=0
prompt crosst2,"Red --> Blue crosstalk"
kill()
load(foldername,filename,last)
threshold_data(donor,acceptor,autod,autob,crosst,crosst2)

end


// Load files only
function load(foldername,name,b)
string foldername
string name	// Name of Files
variable b 	// First and Last Filenumber

variable e,k,c	// Some variables
setdatafolder root:
newdatafolder/s $foldername
make/o/n=1 ChA_all,ChB_all
newpath/Q path1,
pathinfo path1
string path=S_path

string num
// Following converts filenames to  numbers:

for(c=0;c<b;c+=1)
	variable d=c+1
	if(d<2)
		string nametoload=name
	elseif(d<10)
		nametoload=name+"_0"+num2str(d)
	else
		nametoload=name+"_"+num2str(d)
	endif

	// Following to load files

	string namer=num2str(c)
	string loader=path+nametoload
	print loader
	LoadWave/J/D/W/K=0/A loader
	wave wave0,wave1

	//Load all first
	for(e=0;e<=(dimsize(wave0,0));e+=1) // Perform autofluorescence correction

		redimension/N=(k+1) ChA_all,ChB_all

		ChA_all[k]=wave0[e]
		ChB_all[k]=wave1[e]
		k+=1
	
	endfor
	killwaves wave0
	killwaves wave1
endfor

end


// Function to just select data that's above the threshold
function threshold_data(donor,acceptor,autod,autob,crosst,crosst2)
variable donor	// A thresh
variable acceptor // B thresh
variable autod 		// Autofluorescence in channel a
variable autob		// Autofluorescence in channel b
variable crosst		// Crosstalk A --> B
variable crosst2 	// Crosstalk B --> A
wave ChA_all,ChB_all


variable e

variable coinc=0,Arate=0,Brate=0
make/o/N=1 A_events,B_events,B_desynch,a_desynch,zplot,z_desynch				// Make data to store values
for(e=0;e<(dimsize(cha_all,0));e+=1)

	if((ChA_all[e]-autod-crosst2*chb_all[e])>donor && (ChB_all[e]-autob-crosst*cha_all[e])>acceptor)					// select events above thresh - coincident

		redimension/N=(coinc+1) A_events,B_events,zplot
		A_events[coinc]=ChA_all[e]-autod-crosst2*chb_all[e]
		B_events[coinc]=ChB_all[e]-autob-crosst*cha_all[e]
		zplot[coinc]=ln(A_events[coinc]/B_events[coinc])
		
	coinc+=1
	endif
	
	if((ChA_all[e]-autod-crosst2*chb_all[e])>donor)			// Select just A events
	Arate+=1
	endif

	if((ChB_all[e]-autob-crosst*cha_all[e])>acceptor)
	Brate+=1
	endif
endfor

variable d
for(e=0;e<(dimsize(cha_all,0));e+=1)			// For chance coincidence
	variable num2=round(enoise(dimsize(cha_all,0)))
	if(num2<0)
		num2=-1*num2
	endif

	if((ChA_all[e]-autod-crosst2*chb_all[e])>donor && (ChB_all[num2]-autob-crosst*cha_all[num2])>acceptor)

		redimension/N=(d+1) A_desynch,B_desynch,z_desynch
		A_desynch[d]=ChA_all[e]-autod-crosst2*chb_all[e]
		B_desynch[d]=ChB_all[num2]-autob-crosst*cha_all[num2]
		z_desynch[d]=ln(A_desynch[d]/B_desynch[d])

		d+=1
endif
endfor
make/o/n=1 coincident=coinc,donorrate=arate,acceptorrate=brate,desynchrate=d,qvalue
variable q=(coinc-d)/(arate+brate-(coinc-d))
qvalue[0]=q
string q_val=num2str(q)
Make/N=50/O zplot_Hist,zdesynch_hist;DelayUpdate
Histogram/B={-5,0.2,50} zplot,zplot_Hist
Histogram/B={-5,0.2,50} z_desynch,zdesynch_Hist
Display zplot_Hist,zdesynch_hist
Label left "no. of counts";DelayUpdate
Label bottom "ln(Ia/Id)"
ModifyGraph rgb(zplot_Hist)=(0,0,0),rgb(zdesynch_hist)=(21760,21760,21760)
ModifyGraph rgb(zplot_Hist)=(65280,0,0),rgb(zdesynch_hist)=(39168,39168,39168)
Legend/C/N=text0/J/F=0/A=MC "\\s(zplot_Hist) Total events\r\\s(zdesynch_hist) Desynchronised events"
TextBox/C/N=text0/F=0/A=MC "Q = "+q_val+""
ModifyGraph width=283.465,height=198.425,gfSize=20
ModifyGraph mirror=1,minor=1,btLen=3,stLen=2;DelayUpdate
Label bottom "ln(IR/IG)";DelayUpdate
SetAxis/A/N=1 left
ModifyGraph mode=6
CurveFit/NTHR=0 gauss  zplot_Hist /D 
ModifyGraph gFont="Arial"
SetAxis left 0,*

end



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





Macro Preview(first,last)

variable First=1
Prompt First, "First Bin: "			// Set prompt for X param
variable Last=8000
Prompt Last, "Last Bin: "			// Set prompt for X param

preview_disp(first,last)

end

Macro Plot_for_sarah(first,last,bin_time)

variable First=1
Prompt First, "First Bin: "			// Set prompt for X param
variable Last=8000
Prompt Last, "Last Bin: "			// Set prompt for X param
variable bin_time=50e-6
Prompt bin_time, "Bin width (s)"
sarah_disp(first,last,bin_time)

end


function sarah_disp(first,last,bin_time)		// Show bins
variable first,last,bin_time
wave chA_all,chB_all
make/o/n=(dimsize(chA_all,0)) timeax


variable a,b,t=0
for(a=0;a<(dimsize(chA_all,0));a+=1)
timeax[a]=t
t+=bin_time
endfor

first=first*bin_time
last=last*bin_time


display ChB_all vs timeax
SetAxis bottom first,last
ModifyGraph width=425.197,height=283.465,gFont="Arial",gfSize=20
ModifyGraph btLen=5,stLen=4
Label left "Intensity (counts/bin)";DelayUpdate
Label bottom "Time (s)"
ModifyGraph width=850.394


display ChA_all vs timeax
SetAxis bottom first,last
ModifyGraph width=425.197,height=283.465,gFont="Arial",gfSize=20
ModifyGraph btLen=5,stLen=4
ModifyGraph rgb(ChA_all)=(2,39321,1)
Label left "Intensity (counts/bin)";DelayUpdate
Label bottom "Time (s)"
ModifyGraph width=850.394

end



function preview_disp(first,last)		// Show bins
variable first,last
wave chA_all,chB_all

variable a,b

duplicate/o chb_all,chb_all_inv

for(a=0;a<(dimsize(chb_all,0));a+=1)
chb_all_inv[a]=-1*chb_all_inv[a]
endfor

display

AppendToGraph ChA_all
AppendToGraph ChB_all_inv
SetAxis bottom first,last
ModifyGraph width=425.197,height=283.465,gFont="Arial",gfSize=20
ModifyGraph btLen=5,stLen=4
ModifyGraph rgb(ChA_all)=(2,39321,1)
Label left "Counts";DelayUpdate
Label bottom "Bin number"
ModifyGraph width=850.394


end




function crosstalk_test()
wave ChA_all
wave ChB_all
variable thresh=10
variable a,b,c
make/o/n=1 A_ch,B_ch,rat
for(a=0;a<(dimsize(ChA_all,0));a+=1)
if(cha_all[a]>thresh)
	redimension/n=(b+1) A_ch,b_ch,rat
	a_ch[b]=ChA_all[a]
	b_ch[b]=chB_all[a]
	rat[b]=ln(b_ch[b]/a_ch[b])
	
	b+=1

endif
endfor



end




	
	
function maxQ()

variable autod=0.272		
variable autob=1.93		
variable crosst=0.099		
variable crosst2=0 	
wave ChA_all,ChB_all


variable e

variable coinc=0,Arate=0,Brate=0,desynchrate=0

make/o/n=(50,50) matrix
variable a,b,c

for(a=5;a<20;a+=1)
for(b=5;b<20;b+=1)
coinc=0
Arate=0
Brate=0
desynchrate=0

for(e=0;e<(dimsize(cha_all,0));e+=1)

	if((ChA_all[e]-autod-crosst2*chb_all[e])>a && (ChB_all[e]-autob-crosst*cha_all[e])>b)					// select events above thresh - coincident
		
	coinc+=1
	endif
	
	if((ChA_all[e]-autod-crosst2*chb_all[e])>a)			// Select just A events
	Arate+=1
	endif

	if((ChB_all[e]-autob-crosst*cha_all[e])>b)
	Brate+=1
	endif
endfor

variable d
for(e=0;e<(dimsize(cha_all,0));e+=1)			// For chance coincidence
	variable num2=round(enoise(dimsize(cha_all,0)))
	if(num2<0)
		num2=-1*num2
	endif

	if((ChA_all[e]-autod-crosst2*chb_all[e])>a && (ChB_all[num2]-autob-crosst*cha_all[num2])>b)

		
		desynchrate+=1
endif
endfor

variable q=(coinc-desynchrate)/(arate+brate-(coinc-desynchrate))
matrix[a][b]=q

endfor
endfor

wavestats matrix

printf "Maximum Q is %G\r",v_max
printf "Blue threshold is %d\r",v_maxrowloc 
printf "Red threshold is %d\r",v_maxcolloc 
end




macro close_windows()
kill()
end
