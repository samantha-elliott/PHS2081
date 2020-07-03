# PHS2081 - X-ray Powder Diffraction C1
#Student Name: Samantha Elliott
#Student Id: 30119057
# Import relevant functions
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sp
import monashspa.PHS2061 as spa



# Import relevant data, suggested format for import is a .csv
# (comma-seperated values) file. Note that you man need to use
# skip_header to remove column headers in the file.

si_data = np.genfromtxt("Silicon_data.csv", delimiter=",", skip_header=1)
cu_data = np.genfromtxt("Copper_data.csv", delimiter=",", skip_header=1)
al_data = np.genfromtxt("Aluminium_data.csv", delimiter=",", skip_header=1)


# You should unpack any imported data, and create any other required values
TwoTheta = si_data[:,0]
I = si_data[:,1]
#u_2Theta_data = np.zeros_like(2Theta_data)
#u_I_data = np.zeros_like(I_data)

# Create a Plot Si
plt.figure()
plt.errorbar(TwoTheta, I, marker="None", color="black", label="experiment data")
plt.title("Figure: Plot of intensity vs 2θ for Silicon")
plt.xlabel("2θ [degrees]")
plt.ylabel("Intensity")
leg = plt.legend(bbox_to_anchor=(1,1))
plt.savefig('Si_data_graph.png',dpi=600, bbox_extra_artists=(leg,), bbox_inches='tight')
plt.show()

# Find the Peaks Si
# We suggest scipy's find peaks, there are certainly more elegant solutions but this will suffice.
# You will need to provide a "prominence" figure, which is essentially the minimum
# distance above the baseline for something to be a peak
indices, _ = sp.find_peaks(I,prominence=300)
# Next we turn those indices into angles
peaks = np.array([TwoTheta [j] for j in indices])
# and return the full width at half measure, a measure of uncertainty for peaks.
# Note we need to calibrate for the spacing of the x data as the algorithm returns the width in indices
FWHM = (TwoTheta [1]-TwoTheta [0])*np.array(sp.peak_widths(I, indices, rel_height=0.5)[0])
print("Peaks were found at x positions ",peaks," with full widths at half measure of ",FWHM)




#Theta values where peaks are 
Theta = peaks/2
#sin(theta)
SinTheta = np.sin(np.deg2rad(Theta))
#sin^2(theta)
SinThetaSquared = np.square(SinTheta)


SinTheta1 = np.array([0.2464916,0.40178796,0.47101188,0.56755638], dtype=np.float64)

RootN = np.array([1.73,2.83,3.3166,4], dtype=np.float64)

fit_result = spa.linear_fit(SinTheta1, RootN)
y_fit = fit_result.best_fit

plt.figure()
plt.errorbar(SinTheta, RootN, marker="o", color="black", label="experiment data")
plt.title("Figure: Plot of Sin(θ) vs √N for Silicon")
plt.xlabel("Sin(θ) [degrees]")
plt.ylabel("√N")
leg = plt.legend(bbox_to_anchor=(1,1))
plt.savefig('SinTheta vs RootNSi.png',dpi=600, bbox_extra_artists=(leg,), bbox_inches='tight')
plt.show()

#gradient 
fit_parameters=spa.get_fit_parameters(fit_result)
print(fit_parameters)
slope = fit_parameters["slope"]
a = slope*1.54
a = a/2
print("a is", a)


#uncertainty in a 
u_slope = fit_parameters["u_slope"]
u_a = u_slope*1.54
u_a = u_a/2
print("u_a is", u_a)




# COPPER
# You should unpack any imported data, and create any other required values
Theta = cu_data[:,0]
I = cu_data[:,1]
#u_2Theta_data = np.zeros_like(2Theta_data)
#u_I_data = np.zeros_like(I_data)

# Create a Plot Cu
plt.figure()
plt.errorbar(Theta, I, marker="None", color="black", label="experiment data")
plt.title("Figure: Plot of Intensity vs 2θ for Copper")
plt.xlabel("2θ [degrees]")
plt.ylabel("Intensity")
leg = plt.legend(bbox_to_anchor=(1,1))
plt.savefig('Cu_data_graph.png',dpi=600, bbox_extra_artists=(leg,), bbox_inches='tight')
plt.show()

# Find the Peaks Cu
# We suggest scipy's find peaks, there are certainly more elegant solutions but this will suffice.
# You will need to provide a "prominence" figure, which is essentially the minimum
# distance above the baseline for something to be a peak
indices, _ = sp.find_peaks(I,prominence=200)
# Next we turn those indices into angles
peaks = np.array([Theta [j] for j in indices])
# and return the full width at half measure, a measure of uncertainty for peaks.
# Note we need to calibrate for the spacing of the x data as the algorithm returns the width in indices
FWHM = (Theta [1]-Theta [0])*np.array(sp.peak_widths(I, indices, rel_height=0.5)[0])
print("Peaks were found at x positions ",peaks," with full widths at half measure of ",FWHM)


#Theta values where peaks are 
Theta = peaks/2
#sin(theta)
SinTheta = np.sin(np.deg2rad(Theta))
#sin^2(theta)
SinThetaSquared = np.square(SinTheta)

SinTheta = np.array([0.36861133,0.42562136,0.60209376,0.70624236,0.73763097], dtype=np.float64)

RootN = np.array([1.73,2,2.83,3.32,3.46], dtype=np.float64)

fit_result = spa.linear_fit(SinTheta, RootN)
y_fit = fit_result.best_fit

plt.figure()
plt.errorbar(SinTheta, RootN, marker="o", color="black", label="experiment data")
plt.title("Figure: Plot of Sin(θ) vs √N for Copper")
plt.xlabel("Sin(θ) [degrees]")
plt.ylabel("√N")
leg = plt.legend(bbox_to_anchor=(1,1))
plt.savefig('SinTheta vs RootNCo.png',dpi=600, bbox_extra_artists=(leg,), bbox_inches='tight')
plt.show()

#gradient 
fit_parameters=spa.get_fit_parameters(fit_result)
print(fit_parameters)
slope = fit_parameters["slope"]
a = slope*1.54
a = a/2
print("a is", a)


#uncertainty in a 
u_slope = fit_parameters["u_slope"]
u_a = u_slope*1.54
u_a = u_a/2
print("u_a is", u_a)




#Aluminium
# You should unpack any imported data, and create any other required values
Theta = al_data[:,0]
I = al_data[:,1]
#u_2Theta_data = np.zeros_like(2Theta_data)
#u_I_data = np.zeros_like(I_data)

# Create a Plot Cu
plt.figure()
plt.errorbar(Theta, I, marker="None", color="black", label="experiment data")
plt.title("Figure: Plot of rintensity vs 2θ for Aluminium")
plt.xlabel("2θ [degrees]")
plt.ylabel("Intensity")
leg = plt.legend(bbox_to_anchor=(1,1))
plt.savefig('Al_data_graph.png',dpi=600, bbox_extra_artists=(leg,), bbox_inches='tight')
plt.show()

# Find the Peaks Cu
# We suggest scipy's find peaks, there are certainly more elegant solutions but this will suffice.
# You will need to provide a "prominence" figure, which is essentially the minimum
# distance above the baseline for something to be a peak
indices, _ = sp.find_peaks(I,prominence=250)
# Next we turn those indices into angles
peaks = np.array([Theta [j] for j in indices])
# and return the full width at half measure, a measure of uncertainty for peaks.
# Note we need to calibrate for the spacing of the x data as the algorithm returns the width in indices
FWHM = (Theta [1]-Theta [0])*np.array(sp.peak_widths(I, indices, rel_height=0.5)[0])
print("Peaks were found at x positions ",peaks," with full widths at half measure of ",FWHM)

#Theta values where peaks are 
Theta = peaks/2
#sin(theta)
SinTheta = np.sin(np.deg2rad(Theta))
#sin^2(theta)
SinThetaSquared = np.square(SinTheta)

SinTheta = np.array([0.33002017,0.38107038,0.53847668,0.63121744,0.65921458], dtype=np.float64)

RootN = np.array([1.73,2,2.83,3.32,3.46], dtype=np.float64)

fit_result = spa.linear_fit(SinTheta, RootN)
y_fit = fit_result.best_fit

plt.figure()
plt.errorbar(SinTheta, RootN, marker="o", color="black", label="experiment data")
plt.title("Figure: Plot of Sin(θ) vs √N for Aluminium")
plt.xlabel("Sin(θ) [degrees]")
plt.ylabel("√N")
leg = plt.legend(bbox_to_anchor=(1,1))
plt.savefig('SinTheta vs RootNAl.png',dpi=600, bbox_extra_artists=(leg,), bbox_inches='tight')
plt.show()

#gradient 
fit_parameters=spa.get_fit_parameters(fit_result)
print(fit_parameters)
slope = fit_parameters["slope"]
a = slope*1.54
a = a/2
print("a is", a)


#uncertainty in a 
u_slope = fit_parameters["u_slope"]
u_a = u_slope*1.54
u_a = u_a/2
print("u_a is", u_a)