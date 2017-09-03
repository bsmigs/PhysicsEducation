"""
General Numerical Solver for 2D-projectile motion.

author: Brian Smigielski
email: bsmigs@gmail.com
website: http://rfground.wordpress.com
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.integrate as integrate
from scipy import interpolate
from Tkinter import *

class ProjectileMotion:
    
    def __init__(self):
        self.fields = {'mass':1.0, 'g':9.81, 'v0':10, 'theta':45}
        self.unitFields = {'mass':'(kg)', 'g':'(m/s^2)', 'v0':'(m/s)', 'theta':'(deg)'}
        self.dragFields = {'drag coefficient':0.5, 'diameter':0.01, 'air density':1.225} 
        self.unitDragFields = {'drag coefficient':'(number)', 'diameter':'(m)', 'air density':'(kg/m^3)'}
        self.origin = (0, 0)
        self.dt = 0.01
        self.time_elapsed = 0
        self.state = []


    def setValues(self, entries):
        for entry in entries:
            field = entry[0]
            value  = entry[1].get()
            self.fields[field] = value

            try:
                self.fields[field] = float(self.fields[field])
            except ValueError:
                print "Error: You must enter a number for field: ",field.upper()
                abort()

            # make sure all fields > 0. We'll take care of
            # "g" below
            self.fields[field] = np.abs(self.fields[field])

        # make sure "g" is < 0
        self.fields['g'] = -np.abs(self.fields['g'])

        # make sure theta between 0 and 180
        if ((self.fields['theta'] < 0.) | (self.fields['theta'] > 180.)):
            print('Theta must be between 0 and 180 deg')
            abort()
                
        # define initial state
        vels = self.getInitialVelocities()
        self.state = [self.origin[0], vels[0], self.origin[1], vels[1]]

        # set initial conditions
        self.t = np.array([0])
        self.pos = np.array([self.origin[0], self.origin[1]])
        self.v = np.array([vels[0], vels[1]])

        
    def getInitialVelocities(self):
        """compute the current x,y velocities of the projectile"""
        v0x = self.fields['v0']*np.cos(self.fields['theta']*np.pi/180.)
        v0y = self.fields['v0']*np.sin(self.fields['theta']*np.pi/180.)
        return (v0x, v0y)
    

    def statesNoDrag(self):
        """execute one time step of length dt and update state"""
        # x-direction
        self.state[0] += self.state[1]*self.dt
        self.state[1] = self.state[1]
        
        # y-direction
        self.state[2] += self.state[3]*self.dt + 0.5*self.fields['g']*self.dt*self.dt
        self.state[3] += self.fields['g']*self.dt
        
        # advance time
        self.time_elapsed += self.dt

        # cache values
        self.t = np.append(self.t, self.time_elapsed)
        self.pos = np.append(self.pos, [self.state[0], self.state[2]])
        self.v = np.append(self.v, [self.state[1], self.state[3]])

        
    def derivs(self, dummyState, dummyT):
        area = 0.25 * np.pi * self.fields['diameter']
        dragCoeff = 0.5 * self.fields['drag coefficient'] * self.fields['air density'] * area

        # model the quadratic drag force \alpha |v|^2 \hat{v} as \alpha |v| \vec{v}
        # set \beta = \alpha |v|
        dxdt = dummyState[1]
        dydt = dummyState[3]
        
        beta = dragCoeff * np.sqrt(dummyState[1]*dummyState[1] + dummyState[3]*dummyState[3])
        dvxdt = -beta * dummyState[1]
        dvydt = self.fields['g'] - beta * dummyState[3]

        return np.array([dxdt, dvxdt, dydt, dvydt])
        
        
    def computeDerivedQuantities(self):
        # compute momentum
        self.p = self.fields['mass'] * self.v

        # compute kinetic energy
        self.K = 0.5 * self.fields['mass'] * np.sum(np.square(self.v), axis=1)

        # compute potential energy
        self.U = self.fields['mass'] * np.abs(self.fields['g']) * self.pos[:,1]

        # compute force
        self.F = np.diff(self.p, axis=0) / self.dt
        self.F = np.insert(self.F, 0, self.F[0,:], axis=0)


    def evolve(self, usingDragForce):
        t = self.getTimeVec()

        """
        if (usingDragForce == 0):
            for tt in t:
                self.statesNoDrag()
                
            # reshape so that its in rows of x and y components
            self.pos = self.pos.reshape(len(self.pos)/2, 2)
            self.v = self.v.reshape(len(self.v)/2, 2)
        else:
            # integrate to get solutions
            states = integrate.odeint(self.derivs, self.state, t)
            states = np.array(states)
            
            # break out positions/vels from the state vector
            self.pos = states[:,[0,2]]
            self.v = states[:,[1,3]]
            self.t = t
        """

        if (usingDragForce == 0):
            self.fields['diameter'] = 0

        # integrate to get solutions
        states = integrate.odeint(self.derivs, self.state, t)
        states = np.array(states)
            
        # break out positions/vels from the state vector
        self.pos = states[:,[0,2]]
        self.v = states[:,[1,3]]
        self.t = t

        maxRange, maxTime = self.maxRange()
        maxHeight = self.maxHeight()

        self.t = np.arange(0,maxTime,self.dt)
        self.pos = self.pos[0:len(self.t),:]
        self.v = self.v[0:len(self.t),:]

        # compute derived quantities
        self.computeDerivedQuantities()

        return maxTime, maxRange, maxHeight

        

    def maxHeight(self):
        """Get max height of projectile"""
        timeInterp = interpolate.interp1d(self.v[:,1], self.t)
        yposInterp = interpolate.interp1d(self.t, self.pos[:,1])

        peakTime = timeInterp(0.0)
        maxHeight = yposInterp(peakTime)

        return maxHeight
        
        #vels = self.getInitialVelocities()
        #return 0.5*(vels[1]*vels[1])/np.abs(self.fields['g'])

    
    def maxRange(self):
        """Get max range of projectile"""        
        timeInterp = interpolate.interp1d(self.pos[1:len(self.pos),1], self.t[1:len(self.pos)])
        xposInterp = interpolate.interp1d(self.t, self.pos[:,0])

        maxTime = timeInterp(0.0)
        maxRange = xposInterp(maxTime)

        return maxRange, maxTime
        
        #vels = self.getInitialVelocities()
        #return (2.0*vels[0]*vels[1])/np.abs(self.fields['g'])


    def getTimeVec(self):
        # make sure that we round up to the nearest
        # tenth of a second so that we ensure the particle
        # reaches the ground again -- which makes
        # the interpolations well defined
        totalTime = self.totalTime()
        totalTime *= 10
        totalTime = np.ceil(totalTime)
        totalTime /= 10
        
        timeVec = np.arange(0,np.ceil(totalTime),self.dt)
        return timeVec
      
    
    def totalTime(self):
        vels = self.getInitialVelocities()
        return (2.0*vels[1])/np.abs(self.fields['g'])


    def clear(self):
        self.__init__()

##########################################################
""" RUN SIM """
##########################################################

pm = ProjectileMotion()
root = Tk()

row = Frame(root)
label = Label(row, width=15, text='POOP', anchor='w')
label.pack(side=TOP)

entries = []
for key in pm.fields:
    row = Frame(root)
    label = Label(row, width=15, text=key, anchor='w')
    entry = Entry(row)
    entry.insert(END, pm.fields[key])

    row.bind()
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    label.pack(side=LEFT)
    entry.pack(side=LEFT, expand=YES, fill=X)

    # set units
    unitsLabel = Label(row, width=8, text=pm.unitFields[key], anchor='w')
    unitsLabel.pack(side=RIGHT)

    # append values to entries list
    entries.append( (key, entry) )
    
    
# Begin logic to include air resistance
def revealOptions(dragEntries):
    for ent in dragEntries:
        if (CheckVar0.get() == 0):
            ent[1].config(state='disabled')
        else:
            ent[1].config(state='normal')
            
CheckVar0 = IntVar()
C0 = Checkbutton(root, text="Include air resistance", justify=LEFT, variable=CheckVar0, command=lambda: revealOptions(dragEntries))
C0.pack(side=TOP, anchor=W)

dragEntries = []
for key in pm.dragFields:
    row = Frame(root)
    label = Label(row, width=15, text=key, anchor='w')
    ent = Entry(row)
    ent.insert(END, pm.dragFields[key])

    row.bind()
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    label.pack(side=LEFT)
    ent.pack(side=LEFT, expand=YES, fill=X)
    ent.config(state='disabled')

    # set units
    unitsLabel = Label(row, width=8, text=pm.unitDragFields[key], anchor='w')
    unitsLabel.pack(side=RIGHT)

    dragEntries.append((key, ent))
# End logic to include air resistance


CheckVar1 = IntVar()
CheckVar2 = IntVar()
CheckVar3 = IntVar()
CheckVar4 = IntVar()
CheckVar5 = IntVar()
CheckVar6 = IntVar()

C1 = Checkbutton(root, text="Y-position vs. X-position", variable=CheckVar1, justify=LEFT)
C2 = Checkbutton(root, text="X-position vs. Time", variable=CheckVar2, justify=LEFT)
C3 = Checkbutton(root, text="Y-position vs. Time", variable=CheckVar3, justify=LEFT)
C4 = Checkbutton(root, text="Force vs. Time", variable=CheckVar4, justify=LEFT)
C5 = Checkbutton(root, text="Momentum vs. Time", variable=CheckVar5, justify=LEFT)
C6 = Checkbutton(root, text="Energy vs. Time", variable=CheckVar6, justify=LEFT)

C1.pack(side=TOP, anchor=W)
C2.pack(side=TOP, anchor=W)
C3.pack(side=TOP, anchor=W)
C4.pack(side=TOP, anchor=W)
C5.pack(side=TOP, anchor=W)
C6.pack(side=TOP, anchor=W)


def runAnimation():

    # append drag force fields to normal ones
    # no matter what
    allEntries = entries + dragEntries
        
    # read in and set user values
    pm.setValues(allEntries)

    # ask if we are using drag force
    usingDragForce = CheckVar0.get()
    maxTime, maxRange, maxHeight = pm.evolve(usingDragForce)
    
    t_size = pm.t.size
    timeIdxs = np.arange(t_size)

    if (CheckVar1.get() == 1):
        ## set up figure y-position vs. x-position animation
        fig = plt.figure()
        ax = fig.add_subplot(111, autoscale_on=False, \
                            xlim=(0, maxRange+2), ylim=(0, maxHeight+2))
        ax.grid()
        line, = ax.plot(pm.origin[0], pm.origin[1], '-')
        time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
        plt.xlabel("Horizontal position (m)")
        plt.ylabel("Vertical position (m)")
        plt.title("Vertical Position vs. Horizontal Position")

        def y_vs_x(i):
            x = pm.pos[0:i+1, 0]
            y = pm.pos[0:i+1, 1]
            line.set_data(x, y)
            time_text.set_text('time = %.3f (s)' % pm.t[i])
            return line, time_text
    
        ani = animation.FuncAnimation(fig, y_vs_x, frames=timeIdxs,
                                    interval=3000*pm.dt, blit=True, repeat=False)


    ## set up figure x-position vs. time animation
    if (CheckVar2.get() == 1):
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, autoscale_on=False, \
                            xlim=(0, maxTime+0.5), ylim=(0, maxRange+2))
        ax2.grid()
        time_text2 = ax2.text(0.02, 0.95, '', transform=ax2.transAxes)
        line2, = ax2.plot([], [], '-')
        plt.xlabel("Time (s)")
        plt.ylabel("Horizontal position (m)")
        plt.title("Horizontal Position vs. Time")

        def x_vs_t(i):
            t = pm.t[0:i+1]
            x = pm.pos[0:i+1, 0]
            line2.set_data(t, x)
            time_text2.set_text('time = %.3f (s)' % pm.t[i])
            return line2, time_text2

        ani2 = animation.FuncAnimation(fig2, x_vs_t, frames=timeIdxs,
                                    interval=3000*pm.dt, blit=True, repeat=False)


    ## set up figure y-position vs. time animation
    if (CheckVar3.get() == 1):
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111, autoscale_on=False, \
                            xlim=(0, maxTime+0.5), ylim=(0, maxHeight+2))
        ax3.grid()
        time_text3 = ax3.text(0.02, 0.95, '', transform=ax3.transAxes)
        line3, = ax3.plot([], [], '-')
        plt.xlabel("Time (s)")
        plt.ylabel("Vertical position (m)")
        plt.title("Vertical Position vs. Time")

        def y_vs_t(i):
            t = pm.t[0:i+1]
            y = pm.pos[0:i+1, 1]
            line3.set_data(t, y)
            time_text3.set_text('time = %.3f (s)' % pm.t[i])
            return line3, time_text3

        ani3 = animation.FuncAnimation(fig3, y_vs_t, frames=timeIdxs,
                                        interval=3000*pm.dt, blit=True, repeat=False)


    ## set up figure force vs. time animation
    if (CheckVar4.get() == 1):
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

        ani4 = animation.FuncAnimation(fig4, F_vs_t, frames=timeIdxs,
                                        interval=3000*pm.dt, blit=True, repeat=False)



    ## set up figure momentum vs. time animation
    if (CheckVar5.get() == 1):
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

        ani5 = animation.FuncAnimation(fig5, p_vs_t, frames=timeIdxs,
                                        interval=3000*pm.dt, blit=True, repeat=False)


    ## set up figure energy vs. time animation
    if (CheckVar6.get() == 1):
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

        ani6 = animation.FuncAnimation(fig6, E_vs_t, frames=timeIdxs,
                                        interval=3000*pm.dt, blit=True, repeat=False)



    print("Maximum Time = %.3f (s)" % maxTime)
    print("Maximum Range = %.3f (m)" % maxRange)
    print("Maximum Height = %.3f (m)" % maxHeight)
        
    # show all relevant plots and then quit when finished
    # TODO: Clear out pm object so sim can run again
    plt.show()
    #root.quit()
    pm.clear()


# set up buttons to run
b1 = Button(root, text='Run', command=runAnimation)
b1.pack(side=LEFT, padx=5, pady=5)

# set up buttons to quit
b2 = Button(root, text='Quit', command=root.quit)
b2.pack(side=RIGHT, padx=5, pady=5)

# run it
root.mainloop()
