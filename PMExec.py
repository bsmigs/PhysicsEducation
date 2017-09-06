import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import Tkinter as tk
import ProjectileMotion as pm
from WindowConstruction import *
from PMPlots import *

root = tk.Tk()
root.title('Particle 1 Parameter Window')

# Ask eventually if we are looking at 1 or 2 particles
# If 2, we will then compare their performances
# to one another

# create object
pm = pm.ProjectileMotion()

# Add basic fields
entries = loadFields(root, pm.fields, pm.unitFields)
    
# Begin logic to include air resistance
def revealOptions(dragEntries):
    for ent in dragEntries:
        if (CheckVar0.get() == 0):
            ent[1].config(state='disabled')
        else:
            ent[1].config(state='normal')
            
CheckVar0 = tk.IntVar()
C0 = tk.Checkbutton(root, text="Include air resistance", justify=tk.LEFT, variable=CheckVar0, command=lambda: revealOptions(dragEntries))
C0.pack(side=tk.TOP, anchor=tk.W)

dragEntries = loadFields(root, pm.dragFields, pm.unitDragFields, 'disabled')
# End logic to include air resistance

    

typesOfPlots = ['Y-position vs. X-position', 'X-position vs. Time', 'Y-position vs. Time', \
                    'Force vs. Time', 'Momentum vs. Time', 'Energy vs. Time']

checkVarList = []
for plotType in typesOfPlots:
    CheckVar = tk.IntVar()
    checkVarList.append(CheckVar)
    
    C = tk.Checkbutton(root, text=plotType, variable=CheckVar, justify=tk.LEFT)
    C.pack(side=tk.TOP, anchor=tk.W)
    


    
def create_window():
    window = tk.Toplevel(root)

b = tk.Button(root, text="Create new window", command=create_window)
b.pack()




def runAnimation(checkVarList):

    # append drag force fields to normal ones
    # no matter what
    allEntries = entries + dragEntries
    
    # read in and set user values
    pm.setValues(allEntries)


    # ask if we are using drag force
    usingDragForce = CheckVar0.get()
    maxTime, maxRange, maxHeight = pm.evolve(usingDragForce)
    timeIdxs = np.arange(pm.t.size)

    intervalTime = 3000*pm.dt

    
    
    ## set up figure y-position vs. x-position animation
    if (checkVarList[0].get() == 1):
        fig, ax = plt.subplots()
        my_plots = PMPlots(ax, pm)
        my_plots.set_plot_type('y_vs_x')
        ani1 = animation.FuncAnimation(fig, my_plots.set_data, frames=timeIdxs, \
                                    interval=intervalTime, blit=True, repeat=False)

    ## set up figure x-position vs. time animation
    if (checkVarList[1].get() == 1):
        fig, ax = plt.subplots()
        my_plots = PMPlots(ax, pm)
        my_plots.set_plot_type('y_vs_t')
        ani2 = animation.FuncAnimation(fig, my_plots.set_data, frames=timeIdxs, \
                                    interval=intervalTime, blit=True, repeat=False)

    ## set up figure y-position vs. time animation
    if (checkVarList[2].get() == 1):
        fig, ax = plt.subplots()
        my_plots = PMPlots(ax, pm)
        my_plots.set_plot_type('y_vs_t')
        ani3 = animation.FuncAnimation(fig, my_plots.set_data, frames=timeIdxs, \
                                    interval=intervalTime, blit=True, repeat=False)


    ## set up figure force vs. time animation
    if (checkVarList[3].get() == 1):
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(111, autoscale_on=False, \
                            xlim=(0, maxTime+0.5), ylim=(min(pm.F[:,1])-2, max(pm.F[:,0])+2))
        ax4.grid()
        time_text4 = ax4.text(0.02, 0.95, '', transform=ax4.transAxes)
        line4a, = ax4.plot([], [], '-')
        line4b, = ax4.plot([], [], '-')
        plt.xlabel("Time (s)")
        plt.ylabel("Force (N)")
        plt.title("Force vs. Time")
        plt.legend(['Fx (horizontal)', 'Fy (vertical)'])

        def F_vs_t(i):
            t = pm.t[0:i+1]
            Fx = pm.F[0:i+1, 0]
            Fy = pm.F[0:i+1, 1]
            line4a.set_data(t, Fx)
            line4b.set_data(t, Fy)
            time_text4.set_text('time = %.3f (s)' % pm.t[i])
            return line4a, line4b, time_text4

        ani4 = animation.FuncAnimation(fig4, F_vs_t, frames=timeIdxs, \
                                        interval=intervalTime, blit=True, repeat=False)



    ## set up figure momentum vs. time animation
    if (checkVarList[4].get() == 1):
        fig5 = plt.figure()
        ax5 = fig5.add_subplot(111, autoscale_on=False, \
                            xlim=(0, maxTime+0.5), ylim=(min(pm.p[:,1])-2, max(pm.p[:,1])+2))
        ax5.grid()
        time_text5 = ax5.text(0.02, 0.95, '', transform=ax5.transAxes)
        line5a, = ax5.plot([], [], '-')
        line5b, = ax5.plot([], [], '-')
        plt.xlabel("Time (s)")
        plt.ylabel("Momentum (kgm/s)")
        plt.title("Momentum vs. Time")
        plt.legend(['px (horizontal)', 'py (vertical)'])

        def p_vs_t(i):
            t = pm.t[0:i+1]
            px = pm.p[0:i+1, 0]
            py = pm.p[0:i+1, 1]
            line5a.set_data(t, px)
            line5b.set_data(t, py)
            time_text5.set_text('time = %.3f (s)' % pm.t[i])
            return line5a, line5b, time_text5

        ani5 = animation.FuncAnimation(fig5, p_vs_t, frames=timeIdxs, \
                                        interval=intervalTime, blit=True, repeat=False)


    ## set up figure energy vs. time animation
    if (checkVarList[5].get() == 1):
        #E = [pm.K[i]+pm.U[i] for i in xrange(len(pm.K))]
        E = pm.K + pm.U
        
        fig6 = plt.figure()
        ax6 = fig6.add_subplot(111, autoscale_on=False, \
                            xlim=(0, maxTime+0.5), ylim=(min(pm.U)-2, max(E)+2))
        ax6.grid()
        time_text6 = ax6.text(0.02, 0.95, '', transform=ax6.transAxes)
        line6a, = ax6.plot([], [], '-')
        line6b, = ax6.plot([], [], '-')
        line6c, = ax6.plot([], [], '-')
        plt.xlabel("Time (s)")
        plt.ylabel("Energy (J)")
        plt.title("Energy vs. Time")
        plt.legend(['Kinetic', 'Potential', 'Total'])

        def E_vs_t(i):
            t = pm.t[0:i+1]
            K = pm.K[0:i+1]
            U = pm.U[0:i+1]
            Eshort = E[0:i+1]

            line6a.set_data(t, K)
            line6b.set_data(t, U)
            line6c.set_data(t, Eshort)
            time_text6.set_text('time = %.3f (s)' % pm.t[i])
            return line6a, line6b, line6c, time_text6

        ani6 = animation.FuncAnimation(fig6, E_vs_t, frames=timeIdxs, \
                                        interval=intervalTime, blit=True, repeat=False)



    print("Maximum Time = %.3f (s)" % maxTime)
    print("Maximum Range = %.3f (m)" % maxRange)
    print("Maximum Height = %.3f (m)" % maxHeight)
        
    # show all relevant plots and then quit when finished
    plt.show()
    pm.clear()


# set up buttons to run
b1 = tk.Button(root, text='Run', command=lambda l=checkVarList: runAnimation(l))
b1.pack(side=tk.LEFT, padx=5, pady=5)

# set up buttons to quit
b2 = tk.Button(root, text='Quit', command=root.quit)
b2.pack(side=tk.RIGHT, padx=5, pady=5)

# run it
root.mainloop()
