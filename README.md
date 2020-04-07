# Protein-aggregation-analysis
Igor Pro procedures for TCCD and FRET analyses.
Igor Pro by Wavemetrics manual can be found here https://www.wavemetrics.com/products/igorpro/manual
All codes written by Dr. Mathew Horrocks, The University of Edinburgh. Detailed explanation of software use can be found here https://pubs.acs.org/doi/10.1021/acs.analchem.5b01811 in SI. 

As previously by Horrocks \textit{et al.}.\cite{Horrocks2015FastOligomers}
Association quotient (Q) for alignment with DNA and then A{\textbeta} was carried out using custom written Igor Pro script (Wavemetrics) by Dr. Mathew Horrocks. For each experiment, files are loaded into the program and user-input thresholds of 10 counts/bin for both channels to minimise background contribution. Time-bins above these thresholds are saved and photon-bursts values are modified to account for autofluorescence of 0.272 in channel A and 1.93 in B. Cross-talk between A and B $~$10 \% was also accounted for. The modified values are calculated by $I_D = D - A_D$ where I$_D$ is the modified intensity for channel A, D the original intensity and A$_D$ the autofluorescence in that channel. The same is applied to channel B with consideration of the cross-talk by $I_A = A - A_A - C \times D$ where I$_A$ is the modified intensity of channel A, A the original intensity, A$_A$ the autofluorescence and C the cross-talk from A to B. Cross-talk from B to A is negligible. Once these adjustments have been made, the value of Q can be calculated by \textit{Equation 1.5, Section 1.8.1} by consideration of event rates in each channel, the observed rate of coincident events and the estimated rate of coincidence which is calculated by $E = AB\tau$ where {\texttau} is interval time in seconds. The calculated Q-value indicated the fraction of events which are considered coincident. For DNA used in alignment (Section 3.3.3) a minimum Q-value of 0.2 was used.\\
\\
In FRET an Igor Pro script was used. As previously, data for analysis were loaded into the program and user-input thresholds of 10 counts/bin for channel A donor and channel B acceptor implemented as before. Bins with counts greater than the thresholds in both channels, considered oligomeric, were saved for further analysis. Those with counts higher than the threshold in only one channel were saved separately and considered monomeric. As above, consideration was given to autofluorescence and cross-talk between the acceptor and donor channels. For bins deemed to contain oligomeric species, the size of the species was estimated by
\begin{equation*}
    Size = 2 \times \frac{I_D + \frac{1}{\gamma}I_A}{I_{monomer}}
\end{equation*}
where {\textgamma} is the detection efficiency consideration (\textit{Equation 1.3, Section 1.7.2}). These values can also be used to determine the FRET efficiency by \textit{Equation 1.2, Section 1.7.2}. Monomer brightness is estimated by the mean of all non-coincident bursts in the donor channel A. It is also possible to estimate the number of monomers present in each oligomer. This can be done as only donor-labelled monomers can be directly excited which means the total intensity of both donor and acceptor is proportional to the number of donor labels with consideration of adjusted acceptor and donor intensities and the gamma factor. For size estimation, it is assumed that for every acceptor there is a monomer present and so the value is multiplied by 2. To convert this to equivalent aggregate size, the value is then divided by the calculated monomer brightness. This method means that the size is FRET efficiency independent.\\
\\
Fitting of the FRET efficiency and size distributions is also carried out by a custom written Igor Pro procedure. As discussed previously, data are sorted into small, medium-low FRET, medium-high FRET and large size by monomer equivalent calculations. This data is then fitted to a normalised, Gaussian distribution for each species type. Initial estimates for Gaussian fits of medium-mers were set as y-offset, amplitude of peak 1, peak 1 x-centre, peak 1 width, peak 2 amplitude, peak 2 centre, peak 2 width were set as 0, 10, 0.5, 0.1, 10, 0.8, 0.1 respectively. 
