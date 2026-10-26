import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

import lcdatabase as lcdb
import fileretriever
import fitter


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np
from scipy.optimize import curve_fit


#Enter Isotopes portion#
#Fields for each: A, Z (optional, name)
Alist, Zlist = [], []
alabels, zlabels = [], []
isotopelabels = []
isotopes = []
iaea_peaks = None

data = None
bins = 0

data_graph = None
scatter = None
threshold = None
curves = []    #list of every gaussian fit

peak_x_axis = None
peak_y_axis = None

top_3_x = []
totalpeaks = None

#stores peaks for manUAL calibration, corresponding box in GUI.... would be a list of PAIRS of coordinates [(x,y), (x,y)...]
peak_x_axis_man = []
peak_y_axis_man = []

additional_x_fits = []
cid_list = []


#linear regression fields
peak_x_axis_kev = []
lineparams = None
###
class Cursor:
    """
    A cross hair cursor.
    """
    def __init__(self, ax,figure):
        self.ax = ax
        self.fig = figure
        self.vertical_line = ax.axvline(color='k', lw=0.8, ls='--')
        self.vertical_line.set_visible(False)
        # text location in axes coordinates
        self.text = ax.text(0.8, 1.05, '', transform=ax.transAxes)

    def set_cross_hair_visible(self, visible):
        need_redraw = self.vertical_line.get_visible() != visible
        self.vertical_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_redraw

    def on_mouse_move(self, event):
        if thresholdvar.get() == 1:    
            if not event.inaxes:
                need_redraw = self.set_cross_hair_visible(False)
                if need_redraw:
                    #curvecanvas.draw()
                    self.ax.figure.canvas.draw()
            else:
                self.set_cross_hair_visible(True)
                x = event.xdata
                # update the line positions
                self.vertical_line.set_xdata([x])
                #curvecanvas.draw()
                self.ax.figure.canvas.draw()


class Isotope:
    #Initialized by accessing the IAEA livechart.
    def __init__(self, A:int, Z:int,rad_type:str):
        self.name = lcdb.element_name(A,Z)
        self.df = lcdb.df_open(self.name,rad_type)

    #Returns the list of energies sorted by default ‘intensity’
    def energylist(self, param='intensity',paramlist=True):   
        data = self.df[['energy',param]].sort_values(by=param, ascending=False)
        if paramlist:
            return data['energy'].tolist(), data[param].tolist()
        return data['energy'].tolist()



def add_isotopefields(event):   #on <return> key, creates fields/destroys old ones
    #add to column 0 and 1, n rows

    #destroy old data
    if Alist or Zlist:
        for n in range(len(Alist)):
            Alist[n].destroy()
            alabels[n].destroy()
            Zlist[n].destroy()
            zlabels[n].destroy()
        for n in range(len(isotopelabels)):
            isotopelabels[n].destroy()
        Alist.clear()
        Zlist.clear()
        alabels.clear()
        zlabels.clear()
        isotopelabels.clear()

    #default triple alpha source
    A_default = [237,241,244]
    Z_default = [93,95,96]

    
    #generate new fields
    isotopecount = int(countentry.get())
    for n in range(isotopecount):
        Alist.append(tk.Entry(infoframe,width=3))
        Zlist.append(tk.Entry(infoframe,width=3))
        alabels.append(tk.Label(infoframe, text="A:"))
        zlabels.append(tk.Label(infoframe, text="Z:"))
        alabels[n].grid(row=n+1,column=0, padx=1)
        zlabels[n].grid(row=n+1, column=2)
        Alist[n].grid(row=n+1, column=1)
        Zlist[n].grid(row=n+1, column=3)
        Alist[n].insert(0,A_default[n])    #filling default values for A and Z
        Zlist[n].insert(0,Z_default[n])

def confirm():
    global isotopes
    if Alist and Zlist and not isotopes:
        isotoperange = int(countentry.get())
        for n in range(isotoperange):
            A = int(Alist[n].get())
            Z = int(Zlist[n].get())
            #for data manipulation
            isotopes.append(Isotope(A,Z,"a"))

            #for GUI display
            isotopelabels.append(tk.Label(infoframe,text=f'{lcdb.element_name(A,Z)}', padx=10))
            isotopelabels[n].grid(row=n+1,column=5)
        isotopes = sorted(isotopes, key= lambda i: i.energylist(paramlist=False)[0],reverse=True)

def openfile():
    global data
    global bins
    filepath = filedialog.askopenfilename()
    if filepath:
        data, bins = fileretriever.fileopener(filepath)
        data_plot(data,bins)
        filesbutton["state"] = "disabled"

def data_plot(data,bins):
    global data_graph
    x = np.linspace(1,bins, num=bins)
    data_graph = ax.step(x, data, 'k', linewidth=0.5)
    ax.set_xlabel("Bins")
    ax.set_ylabel("Counts")
    curvecanvas.draw()

def scatter_plot():
    global data
    global bins
    global scatter
    global totalpeaks
    global peak_x_axis
    global peak_y_axis
    if data:
        list_peak_x_axis, list_peak_y_axis, totalpeaks = fitter.peakseeker(data,bins,params=[5,10,10])
        peak_x_axis = np.array(list_peak_x_axis)
        peak_y_axis = np.array(list_peak_y_axis)
        if scattervar.get() == 1:
            scatter = ax.scatter(peak_x_axis,peak_y_axis, marker='o', c="r", picker=True)
            fig.canvas.mpl_connect('button_press_event', onclick)
            curvecanvas.draw()
        if scattervar.get() == 0:
            if scatter:
                scatter.remove()
                curvecanvas.draw()

def fitting():
    global data
    global bins
    global curves
    global peak_x_axis
    global peak_y_axis
    global top_3_x
    global isotopes
    global iaea_peaks

    #linear regression globals
    global lineparams
    global peak_x_axis_kev

    if data and isotopes:
        x = np.linspace(1,bins,bins)
        if len(peak_y_axis_man) < len(isotopes):     #if the selected peaks list is shorter than the list of isotopes
            isotopefitdata, peak_locs = fitter.isotopefitter(data, isotopes, x, peak_x_axis, peak_y_axis) #peak_locs refers the x_axis (bin) location
        else:
            isotopefitdata, peak_locs = fitter.isotopefitter(data, isotopes, x, peak_x_axis_man, peak_y_axis_man)
        if curvevar.get() == 1:
            for i in isotopefitdata:
                curves.append(ax.plot(x, fitter.gaussian(x, *i), 'r-'))
                curvecanvas.draw()
            
            #lin fit section
            real_energy_levels = fitter.allpeaksenergylist(isotopes,totalpeaks)
            iaea_peaks = [i.energylist(paramlist=False)[0] for i in isotopes]   #isotopes already sorted
            m,b = fitter.lineareqfinder(peak_locs, iaea_peaks, diagnose=True)
            lineparams = [m,b]
            for i in isotopefitdata:
                top_3_x.append(i[1])

            
            x_axis_kev = fitter.calibrate(x, m, b)
            for i in peak_x_axis:
                peak_x_axis_kev.append(fitter.calibrate(i, m, b))

            #ax2.scatter(peak_x_axis, sorted(real_energy_levels), marker=',', c='r')
            ax2.scatter(peak_x_axis,peak_x_axis_kev, marker=',')   #plots points on line directly. option to enable/disable view?
            ax2.scatter(top_3_x,real_energy_levels[:3], marker=',', c='r')
            ax2.set_xlabel('Bin Number')
            ax2.set_ylabel('Energy (keV)')
            print(peak_x_axis_kev)

            ax2.plot(x, x_axis_kev, 'k')
            ax2.grid()


    if curvevar.get() == 0:
        for i in curves:
            i[0].remove()
            curves.clear()
        curvecanvas.draw()
    add_fit_button["state"] = "active"

def add_fit(data, bins):
    global lineparams
    x = np.linspace(1, bins, bins)
    print(additional_x_fits[0])
    popt, _ = curve_fit(fitter.gaussian, x, data, p0 = [25, additional_x_fits[0],10])
    curves.append(ax.plot(x, fitter.gaussian(x, *popt), 'b-'))
    print(popt)
    print(fitter.calibrate(popt[1], lineparams[0], lineparams[1]))
    curvecanvas.draw()


#event picker for matplotlib points
def onpick(event):
    if curves:
        ind = event.ind
        xtemp = np.array(peak_x_axis_man)
        ytemp = np.array(peak_y_axis_man)
        print('onpick scatter:', ind, xtemp[ind], ytemp[ind])
        additional_x_fits.append(xtemp[ind][0])


def onclick(event):
    global peak_x_axis_man
    global peak_y_axis_man
    global scatter
    global threshold
    #for clearing data before threshold
    global data
    global bins
    global data_graph
    if thresholdvar.get() == 1:

        threshold = event.xdata
        x_temp = []
        y_temp = []
        for i in range(len(peak_x_axis)):
            if peak_x_axis[i] > threshold:
                x_temp.append(float(peak_x_axis[i]))
                y_temp.append(float(peak_y_axis[i]))
        peak_x_axis_man = x_temp
        peak_y_axis_man = y_temp
        scatter.remove()
        scatter = ax.scatter(x_temp, y_temp, marker='o', color='red', picker=True)
        for i in range(0,int(threshold)):
            data[i] = 0
        x = np.linspace(1,bins, num=bins)
        data_graph[0].remove()
        data_graph = ax.step(x, data, 'k', linewidth=0.5)
        curvecanvas.draw()



width, height = 1100, 850

root = tk.Tk()
root.title("input sequence")
root.geometry(f"{width}x{height}")

#isotope information frame, separately stores basic isotope information
infocanvas= tk.Canvas(root,bd=3, height=150,width=175,relief='raised')
infocanvas.grid(row=1,column=0)
infoframe = tk.Frame(infocanvas)
infoframe.place(anchor='nw', relx=0.05, rely=0.05)
namerbutton = tk.Button(root, text="Confirm", command=confirm)
namerbutton.grid(row=2,column=0)

#number of isotopes fields
entryframe = tk.Frame(root)
entryframe.grid(row=0,column=0)
numberprompt = tk.Label(entryframe, text="# of Isotopes:")
countentry = tk.Entry(entryframe, width=2)
countentry.insert(0,3)
numberprompt.grid(row=0, column=0)
countentry.grid(row=0,column=1)
countentry.bind('<Return>', add_isotopefields)

#load file button
filesbutton = tk.Button(root, text="Browse...", command=openfile)
filesbutton.grid(row=3, column=0)

#checkboxes
scattervar = tk.IntVar()
scattertoggle = tk.Checkbutton(root, text="Enable peaks",var = scattervar, onvalue=1, offvalue=0, command=scatter_plot)
scattertoggle.grid(row=4, column=0)

curvevar = tk.IntVar()
curvetoggle = tk.Checkbutton(root, text="Enable fit", var = curvevar, onvalue=1, offvalue=0, command= fitting)
curvetoggle.grid(row=5,column=0)

thresholdvar = tk.IntVar()
thresholdtoggle = tk.Checkbutton(root, text="Threshold Selection", var = thresholdvar, onvalue=1, offvalue=0)
thresholdtoggle.grid(row=6,column=0)

add_fit_button = tk.Button(root, text="Fit Additional Peaks", command= lambda: add_fit(data,bins),state="disabled")
add_fit_button.grid(row=7, column=0)

#graph display
fig = Figure(figsize=(7,7), dpi=100)
line_fig = Figure(figsize=(7,7), dpi=100)
ax = fig.add_subplot()
ax2 = line_fig.add_subplot()

graphs_notebook = ttk.Notebook(root)
graphs_notebook.grid(rowspan=20,row=1,column=1)
curveframe = tk.Frame(graphs_notebook, width=250, height=250)
curveframe.pack(fill="both", expand=True)

lineframe = tk.Frame(graphs_notebook, width=250, height=250)
lineframe.pack(fill="both", expand=True)

curvecanvas = FigureCanvasTkAgg(fig, master=curveframe)
curvecanvas.draw()
curvecanvas.get_tk_widget().pack(side=tk.TOP,fill=tk.BOTH,expand=True)
toolbar = NavigationToolbar2Tk(curvecanvas, root, pack_toolbar=False)
toolbar.update()
toolbar.place(relx=0.35, rely=0.925, anchor='center')

linecanvas = FigureCanvasTkAgg(line_fig, master=lineframe)
linecanvas.draw()
linecanvas.get_tk_widget().pack(side=tk.TOP,fill=tk.BOTH,expand=True)

toolbar = NavigationToolbar2Tk(linecanvas, root, pack_toolbar=False)
toolbar.update()
toolbar.place(relx=0.64, rely=0.925, anchor='center')

graphs_notebook.add(curveframe, text="Spectrum")
graphs_notebook.add(lineframe, text="Fit")

#interactive mpl figure commands
#checks if mouse is moving, if moving and threshold is on, select threshold
cursor = Cursor(ax,fig)
fig.canvas.mpl_connect('motion_notify_event', cursor.on_mouse_move)
fig.canvas.mpl_connect('pick_event', onpick)

root.mainloop()