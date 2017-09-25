import Tkinter as tk
import ProjectileMotion as pm
import matplotlib.pyplot as plt
from PMPlots import *

class ProjectileMotionGUI:

    PLOT_TYPES = ['Y-position vs. X-position', 'X-Position vs. Time', \
                    'Y-Position vs. Time', 'X-Force vs. Time', 'Y-Force vs. Time', \
                    'X-Momentum vs. Time', 'Y-Momentum vs. Time', 'Kinetic Energy vs. Time', \
                    'Potential Energy vs. Time', 'Total Energy vs. Time']

    PLOT_NAME_SHORTCUTS = ['y_vs_x', 'x_vs_t', 'y_vs_t', 'fx_vs_t', 'fy_vs_t', \
                           'px_vs_t', 'py_vs_t', 'K_vs_t', 'U_vs_t', 'E_vs_t']

    
    def __init__(self, master):
        self.master = master
        master.title("Projectile Motion Dynamics") # set title for main window

        ### BASIC PARAMETER FRAME ###
        label = tk.Label(text="Basic Parameters", fg="blue", width=15, anchor='w')
        label.pack(side=tk.TOP)
        label.config(font=("Arial", 18))

        # projectile motion object
        self.pm = pm.ProjectileMotion()

        # Populate basic parameter fields for projectile motiom
        self.basicParamEntries = self.loadFields(self.pm.basicParams, self.pm.basicParamsUnits)

        ### AIR RESISTANCE FRAME ###
        label = tk.Label(text="Drag Parameters", fg="blue", width=15, anchor='w')
        label.pack(side=tk.TOP)
        label.config(font=("Arial", 18))

        # Button to include air resistance (or not)
        self.var = tk.IntVar()
        cb = tk.Checkbutton(text="Include air resistance", justify=tk.LEFT, variable=self.var, command=self.revealOptions)
        cb.pack(side=tk.TOP, anchor=tk.W)

        # Populate drag/air-resistance fields with a "disabled" default
        self.dragParamEntries = self.loadFields(self.pm.dragParams, self.pm.dragParamsUnits, 'disabled')

        label = tk.Label(text="Generate Plots", fg="blue", width=15, anchor='w')
        label.pack(side=tk.TOP)
        label.config(font=("Arial", 18))
        
        # Add plot type check buttons
        self.plotChoices()
        
        # run simulation button (the command should load final choices and then run sim)
        rb = tk.Button(master, text='Run Simulation', command=self.runSimulation)
        rb.pack(side=tk.LEFT, padx=5, pady=5)

        # quit buttons
        qb = tk.Button(master, text='Quit Simulation', command=master.quit)
        qb.pack(side=tk.RIGHT, padx=5, pady=5)


    def loadFields(self, fields, unitFields, currState='normal'):
        #frame.bind()
        
        entries = []
        for key in fields:
            # define new frame (row) for each field
            frame = tk.Frame(self.master)
            frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
            
            # label for the current parameter
            label = tk.Label(frame, width=15, text=key, anchor='w')
            label.pack(side=tk.LEFT)

            # entry object which ingests a value for current parameter
            entry = tk.Entry(frame)
            entry.insert(tk.END, fields[key])
            entry.pack(side=tk.LEFT, expand=tk.YES, fill=tk.X)
            entry.config(state=currState)

            # set units
            unitsLabel = tk.Label(frame, width=8, text=unitFields[key], anchor='w')
            unitsLabel.pack(side=tk.RIGHT)

            # append values to entries list
            entries.append( (key, entry) )

        return entries

            
    def revealOptions(self):
        for ent in self.dragParamEntries:
            if (self.var.get() == 0):
                ent[1].config(state='disabled')
            else:
                ent[1].config(state='normal')


    def plotChoices(self):
        # ordering of self.checkVarList same as PLOT_OPTIONS
        self.checkVarList = []        
        for plotType in self.PLOT_TYPES:
            CheckVar = tk.IntVar()
            self.checkVarList.append(CheckVar)
            
            cb = tk.Checkbutton(text=plotType, variable=CheckVar, justify=tk.LEFT)
            cb.pack(side=tk.TOP, anchor=tk.W)


    def runSimulation(self):
        # This method will take last recorded values
        # and check box parameters and pass them
        # to the projectile motion object in order
        # for it to do computations. Then it will make plots.
        # Essentially, it is the workhorse method

        # append drag force fields to normal ones
        #allEntries = self.basicParamEntries + self.dragParamEntries

        # read in and set user values
        self.pm.setValues(self.basicParamEntries, 'basic')
        self.pm.setValues(self.dragParamEntries, 'drag')

        # ask if we are using drag force
        maxTime, maxRange, maxHeight = self.pm.evolve(self.var.get())
        print("Maximum Time = %.3f (s)" % maxTime)
        print("Maximum Range = %.3f (m)" % maxRange)
        print("Maximum Height = %.3f (m)" % maxHeight)

        ani_list = []
        for kk in np.arange(0, len(self.checkVarList)):
            if (self.checkVarList[kk].get() == 1):
                my_plots = PMPlots(self.pm, self.PLOT_NAME_SHORTCUTS[kk])
                ani_list.append(my_plots.runSimulation())

        # show all relevant plots and then clear params when finished
        plt.show()
        self.pm.clear()


        
### MAIN EXECUTABLE ###
if __name__ == "__main__":
    root = tk.Tk()
    my_gui = ProjectileMotionGUI(root)
    root.mainloop()
