from packages import *

#init master graph
x_axis = np.linspace(1,bins,num=bins)

def gaussian(x, a, b, c):
    return a*np.exp( -((x-b)**2)/(2*c**2) )

def peakseeker(data,bins, params=[5,10,10]):
#Detects major peaks in raw data graph. Returns largest peaks up to int "major_peaks_amount",
#and total number of peaks.
    #init peak detection variables
    peak_x_axis = []
    peak_y_axis = []

    multiplier = int(bins/2048)   #originally made peakseeker with 2048 in mind. Change in bin size therefore scales params.



    peaks, _ = scipy.signal.find_peaks(data, height=params[0], prominence=params[1]*multiplier, distance=params[2]*multiplier)
    for peak in peaks:
        peak_x_axis.append(int(peak))    #appends bin location of peak to x-axis (bins)
        peak_y_axis.append(data[peak])   #appends corresponding peak value to y-axis (values)

    return peak_x_axis, peak_y_axis, len(peak_x_axis)

def peaksorter(x_axis, y_axis):
    sortedlist = list(sorted(zip(x_axis, y_axis), key=lambda i: i[1], reverse=True))
    lpeaks_x, lpeaks_y = list(zip(*sortedlist))
    lpeaks_x, lpeaks_y = np.asarray(lpeaks_x), np.asarray(lpeaks_y)
    return lpeaks_x, lpeaks_y

def energyaverage(energies, intensities, threshold:list=[4,25,54.9],iter=2,diagnose=True):
#Averages peaks that are within a predefined threshold.
#Intended to average energy/intensity peaks. “threshold” is a list containing 2 criteria:
#"5" is edif threshold regardless of idif, 25 is edif and 54.9 is idif in %.
        minthresh, difthresh, idifthresh = threshold[0],threshold[1],threshold[2]
        energy_list = list(zip(energies,intensities))
        setlist = []
        avg_energies, avg_intensities = energies, intensities
        for i in range(len(energy_list)):
            for j in range(len(energy_list)):
                dif = abs(energy_list[i][0]-energy_list[j][0])
                idif = abs(((energy_list[i][1]-energy_list[j][1])/energy_list[i][1])*100)
                if dif < minthresh and dif !=0 or (dif < difthresh and dif !=0 and idif < idifthresh):
                #If edif is less than 5keV; or is less than 25 AND also has idif of 54.9
                    a = set({i,j})
                    if a not in setlist:
                        setlist.append(a)
                        p1, p2 = energy_list[i][1], energy_list[j][1]
                        pk1, pk2 = energy_list[i][0], energy_list[j][0]
                        if diagnose and (p1+p2 > 10):
                            print(f"    i% percent difference:{idif: .1f}%")
                            print(f"     Index {i}            Index {j}")
                            print(f"    {energy_list[i]}, {energy_list[j]} >>> {(p1/(p1+p2) * pk1 + p2/(p1+p2) * pk2): .1f}, {(p1+p2): .3f}%\n")
                        avg_energy = round(p1/(p1+p2) * pk1 + p2/(p1+p2) * pk2, 1)
                        avg_energies.append(avg_energy)
                        avg_prob = round((p1+p2), 3)
                        avg_intensities.append(avg_prob)
        if diagnose:
            print("Iteration over")
        flattenedlist = list(set([x for pair in setlist for x in pair]))
        for i in flattenedlist:
            avg_energies.remove(energy_list[i][0])
            avg_intensities.remove(energy_list[i][1])
        avg_energylist = list(sorted(zip(energies,intensities), key=lambda i: i[1], reverse=True))
        return list(zip(*avg_energylist))

###do not generalize to all isotopes, specific to neptunium
def iterative_average(energies:list, intensities:list, threshold=[5,25,54.9], iter = 2):
    if iter == 1:
        return energyaverage(energies, intensities, threshold)
    else:
        avg_energies, avg_intensities = energyaverage(energies, intensities, threshold)
        return iterative_average(list(avg_energies), list(avg_intensities), threshold, iter = iter-1)

def isotopefitter(data,isotopes:list, x_axis, peak_x_axis ,peak_y_axis ,diagnose=False):
    #Uses curve_fit on each isotope peak to get optimal gaussian fits.
    #Takes peak coordinates and sorts them.
    #Returns popt for each isotope (list) and bin location of major peak in each isotope (list)
    peak_locs = []
    poptarray = []

    major_peaks_amount = major_peaks_input #input variable
    largest_peaks = list(sorted(zip(peak_x_axis, peak_y_axis),
                                key=lambda i: i[1], reverse=True)[:major_peaks_amount])
    #largest_peaks creates list of the ordered pair (peak_x_axis, peak_y_axis),
    #then sorts it by largest y-axis using lambda to point to index=1 and reverse to sort from highest
    
    lpeaks_x, lpeaks_y = list(zip(*largest_peaks))      #init largest peaks coordinates
    count = 0
    for i in isotopes:
        name = i.name
        popt, pcov = curve_fit(gaussian, x_axis, data, p0 = [25, lpeaks_x[count],10])
        if diagnose:
            print(f"{name}:")
            print(f"Bin Index: {popt[1]} +- {popt[2]}")
        peak_locs.append(popt[1])
        poptarray.append(popt)
        count+=1
    if diagnose:
        print(f"Peak coordinates: ({lpeaks_x} , {lpeaks_y})")
    print("\n")
    return poptarray, peak_locs

def lineareqfinder(peak_locations, top_peaks,diagnose=False):
    peak_locs = sorted(peak_locations)
    top_peaks = sorted(top_peaks)
    energy_peaks = []
    linear_fit = linregress(peak_locs, top_peaks)
    m, b = linear_fit.slope, linear_fit.intercept
    linparams = [m,b]
    if diagnose:
        print(f"y = {m}x + {b}")
    return linparams

def calibrate(axis, m, b):
    #transforms axis based on m and b
    return m*axis + b

def allpeaksenergylist(Isotopes:list, detectedpeaks:int, by="intensity",diagnose=False):
    #Takes energy list of all given isotopes, and takes all detected peaks detected by peakseeker.
    #Default energy list is sorted by “intensity”.
    #Sorts entire list up to "detectedpeaks" (peaks detected in graph)
    #Intended to be used to compare the linear plot to IAEA data.
    xlist = []
    ylist = []
    y = by
    for i in Isotopes:
        energy, y = i.energylist(param=by, paramlist=True)
        energy, y = iterative_average(energy, y, iter=2)
        xlist += energy
        ylist += y
        ylist = [0 if np.isnan(i) else i for i in ylist]
    if diagnose:
        print(f"XY list: {list(zip(xlist,ylist))}")
    sortedenergies = sorted(list(zip(xlist,ylist)), key=lambda i: i[1], reverse=True)[:detectedpeaks]
    energies, intensities = zip(*sortedenergies)
    if diagnose:
        print(f"energy v. intensity: {list(sortedenergies)}")
    return energies

