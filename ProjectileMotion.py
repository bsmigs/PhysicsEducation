"""
General Numerical Solver for 2D-projectile motion.

author: Brian Smigielski
email: bsmigs@gmail.com
website: http://rfground.wordpress.com
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from Tkinter import *

class ProjectileMotion:
    """Projectile Motion Class

    init_params is [x0, y0, v0, theta], where
    x0 = initial horizontal position (m)
    y0 = initial vertical position (m)
    v0 = initial speed (m/s)
    theta = initial take off angle (deg)

    """
    def __init__(self):
        self.fields = {'mass':1.0, 'g':9.81, 'v0':10, 'theta':45}
        self.origin = (0, 0)
        self.entries = []
        self.dt = 0.01
        self.time_elapsed = 0
        self.state = []

        # set up relevant arrays
        self.x = [self.origin[0]] # horizontal position (m)
        self.y = [self.origin[1]] # vertical position (m)
        self.vx = [] # horizontal speed (m/s)
        self.vy = [] # vertical speed (m/s)
        self.t = [0] # time (s)
        self.Fx = [0] # horizontal force (N)
        self.Fy = [0] # vertical force (N)
        self.px = [] # horizontal momentum (kgm/s)
        self.py = [] # vertical momentum (kgm/s)
        self.K = [] # kinetic energy (J)
        self.U = [] # potential energy (J)

        
    def getInitialVelocities(self):
        """compute the current x,y velocities of the projectile"""
        v0x = self.fields['v0']*np.cos(self.fields['theta']*np.pi/180.)
        v0y = self.fields['v0']*np.sin(self.fields['theta']*np.pi/180.)
        return (v0x, v0y)


    def cacheValues(self):
        # cache values
        self.x.append(self.state[0])
        self.y.append(self.state[2])
        self.vx.append(self.state[1])
        self.vy.append(self.state[3])
        self.t.append(self.time_elapsed)

        ##############################
        # Compute derived quantities
        ##############################
        # force
        Fx = self.fields['mass'] * (self.vx[-1] - self.vx[-2]) / self.dt
        Fy = self.fields['mass'] * (self.vy[-1] - self.vy[-2]) / self.dt
        self.Fx.append(Fx)
        self.Fy.append(Fy)

        # momentum
        px = self.fields['mass'] * self.vx[-1]
        py = self.fields['mass'] * self.vy[-1]
        self.px.append(px)
        self.py.append(py)

        # energy
        K = 0.5 * self.fields['mass'] * (np.dot(self.vx[-1], self.vx[-1]) + np.dot(self.vy[-1], self.vy[-1]))
        U = self.fields['mass'] * np.abs(self.fields['g']) * self.y[-1]
        self.K.append(K)
        self.U.append(U)
    

    def step(self):
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
        self.cacheValues()
        

    def evolve(self):
        t = self.getTimeVec()
        for tt in t:
            self.step()
        

    def setValues(self):
        for entry in self.entries:
            field = entry[0]
            value  = entry[1].get()
            self.fields[field] = value

            try:
                self.fields[field] = float(self.fields[field])
                if (field == "g"):
                    # make sure "g" is < 0
                    self.fields[field] = -np.abs(self.fields[field])
                
            except ValueError:
                print "Error: You must enter a number for field: ",field.upper()
                abort()

        # define initial state
        vels = self.getInitialVelocities()
        self.vx.append(vels[0])
        self.vy.append(vels[1])
        self.px.append(self.fields['mass']*vels[0])
        self.py.append(self.fields['mass']*vels[1])
        self.K.append(0.5*self.fields['mass']*(vels[0]*vels[0] + vels[1]*vels[1]))
        self.U.append(self.fields['mass']*np.abs(self.fields['g'])*self.y[0])
        self.state = [self.origin[0], vels[0], self.origin[1], vels[1]]
    

    def maxHeight(self):
        """Get max height of projectile"""
        vels = self.getInitialVelocities()
        return 0.5*(vels[1]*vels[1])/np.abs(self.fields['g'])

    
    def maxRange(self):
        """Get max range of projectile"""
        vels = self.getInitialVelocities()
        return (2.0*vels[0]*vels[1])/np.abs(self.fields['g'])


    def getTimeVec(self):
        timeVec = np.arange(0,self.totalTime(),self.dt)
        return timeVec
      
    
    def totalTime(self):
        vels = self.getInitialVelocities()
        return (2.0*vels[1])/np.abs(self.fields['g'])

##########################################################
""" RUN SIM """
##########################################################

pm = ProjectileMotion()
root = Tk()

CheckVar0 = IntVar()

for key in pm.fields:
    row = Frame(root)
    label = Label(row, width=15, text=key, anchor='w')
    entry = Entry(row)
    entry.insert(END, pm.fields[key])

    row.bind()
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    label.pack(side=LEFT)
    entry.pack(side=RIGHT, expand=YES, fill=X)

    (pm.entries).append((key, entry))

CheckVar1 = IntVar()
CheckVar2 = IntVar()
CheckVar3 = IntVar()
CheckVar4 = IntVar()
CheckVar5 = IntVar()
CheckVar6 = IntVar()

C1 = Checkbutton(root, text = "Y-position vs. X-position", variable = CheckVar1, justify=LEFT)
C2 = Checkbutton(root, text = "X-position vs. Time", variable = CheckVar2, justify=LEFT)
C3 = Checkbutton(root, text = "Y-position vs. Time", variable = CheckVar3, justify=LEFT)
C4 = Checkbutton(root, text = "Force vs. Time", variable = CheckVar4, justify=LEFT)
C5 = Checkbutton(root, text = "Momentum vs. Time", variable = CheckVar5, justify=LEFT)
C6 = Checkbutton(root, text = "Energy vs. Time", variable = CheckVar6, justify=LEFT)

C1.pack(side=TOP, anchor=W)
C2.pack(side=TOP, anchor=W)
C3.pack(side=TOP, anchor=W)
C4.pack(side=TOP, anchor=W)
C5.pack(side=TOP, anchor=W)
C6.pack(side=TOP, anchor=W)


def runAnimation():
    # read in user values
    pm.setValues()
    # evolve the states for all time
    pm.evolve()
    t_size = pm.getTimeVec().size

    if (CheckVar1.get() == 1):
        ## set up figure y-position vs. x-position animation
        fig = plt.figure()
        ax = fig.add_subplot(111, autoscale_on=False, \
                            xlim=(0, pm.maxRange()+2), ylim=(0, pm.maxHeight()+2))
        ax.grid()
        line, = ax.plot(pm.origin[0], pm.origin[1], '-')
        time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
        plt.xlabel("Horizontal position (m)")
        plt.ylabel("Vertical position (m)")
        plt.title("Vertical Position vs. Horizontal Position")

        def y_vs_x(i):
            x = pm.x[0:i+1]
            y = pm.y[0:i+1]
            line.set_data(x, y)
            time_text.set_text('time = %.3f (s)' % pm.t[i])
            return line, time_text
    
        ani = animation.FuncAnimation(fig, y_vs_x, frames=np.arange(t_size),
                                    interval=3000*pm.dt, blit=True, repeat=False)


    ## set up figure x-position vs. time animation
    if (CheckVar2.get() == 1):
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, autoscale_on=False, \
                            xlim=(0, pm.totalTime()+0.5), ylim=(0, pm.maxRange()+2))
        ax2.grid()
        time_text2 = ax2.text(0.02, 0.95, '', transform=ax2.transAxes)
        line2, = ax2.plot([], [], '-')
        plt.xlabel("Time (s)")
        plt.ylabel("Horizontal position (m)")
        plt.title("Horizontal Position vs. Time")

        def x_vs_t(i):
            t = pm.t[0:i+1]
            x = pm.x[0:i+1]
            line2.set_data(t, x)
            time_text2.set_text('time = %.3f (s)' % pm.t[i])
            return line2, time_text2

        ani2 = animation.FuncAnimation(fig2, x_vs_t, frames=np.arange(t_size),
                                    interval=3000*pm.dt, blit=True, repeat=False)


    ## set up figure y-position vs. time animation
    if (CheckVar3.get() == 1):
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111, autoscale_on=False, \
                            xlim=(0, pm.totalTime()+0.5), ylim=(0, pm.maxHeight()+2))
        ax3.grid()
        time_text3 = ax3.text(0.02, 0.95, '', transform=ax3.transAxes)
        line3, = ax3.plot([], [], '-')
        plt.xlabel("Time (s)")
        plt.ylabel("Vertical position (m)")
        plt.title("Vertical Position vs. Time")

        def y_vs_t(i):
            t = pm.t[0:i+1]
            y = pm.y[0:i+1]
            line3.set_data(t, y)
            time_text3.set_text('time = %.3f (s)' % pm.t[i])
            return line3, time_text3

        ani3 = animation.FuncAnimation(fig3, y_vs_t, frames=np.arange(t_size),
                                        interval=3000*pm.dt, blit=True, repeat=False)


    ## set up figure force vs. time animation
    if (CheckVar4.get() == 1):
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(111, autoscale_on=False, \
                            xlim=(0, pm.totalTime()+0.5), ylim=(min(pm.Fy)-2, max(pm.Fy)+2))
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
            Fx = pm.Fx[0:i+1]
            Fy = pm.Fy[0:i+1]
            line4a.set_data(t, Fx)
            line4b.set_data(t, Fy)
            time_text4.set_text('time = %.3f (s)' % pm.t[i])
            return line4a, line4b, time_text4

        ani4 = animation.FuncAnimation(fig4, F_vs_t, frames=np.arange(t_size),
                                        interval=3000*pm.dt, blit=True, repeat=False)



    ## set up figure force vs. time animation
    if (CheckVar5.get() == 1):
        fig5 = plt.figure()
        ax5 = fig5.add_subplot(111, autoscale_on=False, \
                            xlim=(0, pm.totalTime()+0.5), ylim=(min(pm.py)-2, max(pm.py)+2))
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
            px = pm.px[0:i+1]
            py = pm.py[0:i+1]
            line5a.set_data(t, px)
            line5b.set_data(t, py)
            time_text5.set_text('time = %.3f (s)' % pm.t[i])
            return line5a, line5b, time_text5

        ani5 = animation.FuncAnimation(fig5, p_vs_t, frames=np.arange(t_size),
                                        interval=3000*pm.dt, blit=True, repeat=False)


    ## set up figure energy vs. time animation
    if (CheckVar6.get() == 1):
        E = [pm.K[i]+pm.U[i] for i in xrange(len(pm.K))]
        
        fig6 = plt.figure()
        ax6 = fig6.add_subplot(111, autoscale_on=False, \
                            xlim=(0, pm.totalTime()+0.5), ylim=(min(pm.U)-2, max(E)+2))
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

        ani6 = animation.FuncAnimation(fig6, E_vs_t, frames=np.arange(t_size),
                                        interval=3000*pm.dt, blit=True, repeat=False)



    print("Maximum Time = %.3f (s)" % pm.totalTime())
    print("Maximum Range = %.3f (m)" % pm.maxRange())
    print("Maximum Height = %.3f (m)" % pm.maxHeight())
        
    # show all relevant plots and then
    # quit when finished
    # TODO: Allow root to "reset" so sim
    # can be run again
    plt.show()
    root.quit()


# set up buttons to run
b1 = Button(root, text='Run', command=runAnimation)
b1.pack(side=LEFT, padx=5, pady=5)

# set up buttons to quit
b2 = Button(root, text='Quit', command=root.quit)
b2.pack(side=RIGHT, padx=5, pady=5)

# run it
root.mainloop()




