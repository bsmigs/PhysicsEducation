import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class PMPlots:

    SCALE_FACTOR = 1.05

    def __init__(self, pm, plot_type):

        self.pm = pm
        self.intervalTime = 3000*pm.dt
        self.plot_type = plot_type
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot([], [], 'k-')
        self.time_text = self.ax.text(0.02, 0.95, '', transform=self.ax.transAxes)
        self.ax.grid(True)


    def set_data(self, i, x, y):
        t = self.pm.t
        
        self.line.set_data(x[0:i+1], y[0:i+1])
        self.time_text.set_text('time = %.3f (s)' % t[i])
        return self.line, self.time_text


    def runSimulation(self):

        t = self.pm.t
        timeIdxs = np.arange(t.size)
        
        if (self.plot_type == 'y_vs_x'):
            self.ax.set_xlabel("Horizontal position (m)")
            self.ax.set_ylabel("Vertical position (m)")
            self.ax.set_title("Vertical Position vs. Horizontal Position")

            x = self.pm.pos[:,0]
            y = self.pm.pos[:,1]

            min_x = np.min(x)
            min_y = np.min(y)

            print min_x,min_y

            max_x = np.max(x)
            max_y = np.max(y)

            min_x, max_x = self.findMinMax(min_x, max_x)
            min_y, max_y = self.findMinMax(min_y, max_y)

            self.ax.set_xlim(min_x, max_x)
            self.ax.set_ylim(min_y, max_y)

            ani = animation.FuncAnimation(self.fig, self.set_data, frames=timeIdxs, \
                                              fargs=(x,y), interval=self.intervalTime, \
                                              blit=True, repeat=False)

        elif (self.plot_type == 'x_vs_t'):
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Horizontal position (m)")
            self.ax.set_title("Horizontal Position vs. Time")

            t = self.pm.t
            x = self.pm.pos[:,0]
            y = self.pm.pos[:,1]

            min_pos = np.min([np.min(x), np.min(y)])
            max_pos = np.max([np.max(x), np.max(y)])
            min_pos, max_pos = self.findMinMax(min_pos, max_pos)
            
            self.ax.set_xlim(0, np.max(t)*self.SCALE_FACTOR)
            self.ax.set_ylim(min_pos, max_pos)

            ani = animation.FuncAnimation(self.fig, self.set_data, frames=timeIdxs, \
                                              fargs=(t,x), interval=self.intervalTime, \
                                              blit=True, repeat=False)

        elif (self.plot_type == 'y_vs_t'):
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Vertical position (m)")
            self.ax.set_title("Vertical Position vs. Time")

            t = self.pm.t
            x = self.pm.pos[:,0]
            y = self.pm.pos[:,1]

            min_pos = np.min([np.min(x), np.min(y)])
            max_pos = np.max([np.max(x), np.max(y)])
            min_pos, max_pos = self.findMinMax(min_pos, max_pos)
            
            self.ax.set_xlim(0, np.max(t)*self.SCALE_FACTOR)
            self.ax.set_ylim(0, max_pos)

            ani = animation.FuncAnimation(self.fig, self.set_data, frames=timeIdxs, \
                                              fargs=(t,y), interval=self.intervalTime, \
                                              blit=True, repeat=False)

        elif (self.plot_type == 'fx_vs_t'):
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Horizontal force (N)")
            self.ax.set_title("Horizontal Force vs. Time")

            t = self.pm.t
            fx = self.pm.F[:,0]
            fy = self.pm.F[:,1]

            min_force = np.min([np.min(fx), np.min(fy)])
            max_force = np.max([np.max(fx), np.max(fy)])
            min_force, max_force = self.findMinMax(min_force, max_force)
            
            self.ax.set_xlim(0, np.max(t)*self.SCALE_FACTOR)
            self.ax.set_ylim(min_force, max_force)

            ani = animation.FuncAnimation(self.fig, self.set_data, frames=timeIdxs, \
                                              fargs=(t,fx), interval=self.intervalTime, \
                                              blit=True, repeat=False)

        elif (self.plot_type == 'fy_vs_t'):
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Vertical force (N)")
            self.ax.set_title("Vertical Force vs. Time")

            t = self.pm.t
            fx = self.pm.F[:,0]
            fy = self.pm.F[:,1]

            min_force = np.min([np.min(fx), np.min(fy)])
            max_force = np.max([np.max(fx), np.max(fy)])
            min_force, max_force = self.findMinMax(min_force, max_force)
            
            self.ax.set_xlim(0, np.max(t)*self.SCALE_FACTOR)
            self.ax.set_ylim(min_force, max_force)

            ani = animation.FuncAnimation(self.fig, self.set_data, frames=timeIdxs, \
                                              fargs=(t,fy), interval=self.intervalTime, \
                                              blit=True, repeat=False)

        elif (self.plot_type == 'px_vs_t'):
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Horizontal momentum (kgm/s)")
            self.ax.set_title("Horizontal Momentum vs. Time")

            t = self.pm.t
            px = self.pm.p[:,0]
            py = self.pm.p[:,1]

            min_momentum = np.min([np.min(px), np.min(py)])
            max_momentum = np.max([np.max(px), np.max(py)])
            min_momentum, max_momentum = self.findMinMax(min_momentum, max_momentum)
            
            self.ax.set_xlim(0, np.max(t)*self.SCALE_FACTOR)
            self.ax.set_ylim(min_momentum, max_momentum)

            ani = animation.FuncAnimation(self.fig, self.set_data, frames=timeIdxs, \
                                              fargs=(t,px), interval=self.intervalTime, \
                                              blit=True, repeat=False)

        elif (self.plot_type == 'py_vs_t'):
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Vertical momentum (N)")
            self.ax.set_title("Vertical Momentum vs. Time")

            t = self.pm.t
            px = self.pm.p[:,0]
            py = self.pm.p[:,1]

            min_momentum = np.min([np.min(px), np.min(py)])
            max_momentum = np.max([np.max(px), np.max(py)])
            min_momentum, max_momentum = self.findMinMax(min_momentum, max_momentum)
            
            self.ax.set_xlim(0, np.max(t)*self.SCALE_FACTOR)
            self.ax.set_ylim(min_momentum, max_momentum)

            ani = animation.FuncAnimation(self.fig, self.set_data, frames=timeIdxs, \
                                              fargs=(t,py), interval=self.intervalTime, \
                                              blit=True, repeat=False)

        elif (self.plot_type == 'K_vs_t'):
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Kinetic energy (J)")
            self.ax.set_title("Kinetic Energy vs. Time")

            t = self.pm.t
            K = self.pm.K
            U = self.pm.U
            E = K+U

            min_energy = np.min([np.min(K), np.min(U), np.min(E)])
            max_energy = np.max([np.max(K), np.max(U), np.max(E)])
            min_energy, max_energy = self.findMinMax(min_energy, max_energy)
            
            self.ax.set_xlim(0, np.max(t)*self.SCALE_FACTOR)
            self.ax.set_ylim(min_energy, max_energy)

            ani = animation.FuncAnimation(self.fig, self.set_data, frames=timeIdxs, \
                                              fargs=(t,K), interval=self.intervalTime, \
                                              blit=True, repeat=False)

        elif (self.plot_type == 'U_vs_t'):
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Potential energy (J)")
            self.ax.set_title("Potential Energy vs. Time")

            t = self.pm.t
            K = self.pm.K
            U = self.pm.U
            E = K+U

            min_energy = np.min([np.min(K), np.min(U), np.min(E)])
            max_energy = np.max([np.max(K), np.max(U), np.max(E)])
            min_energy, max_energy = self.findMinMax(min_energy, max_energy)
            
            self.ax.set_xlim(0, np.max(t)*self.SCALE_FACTOR)
            self.ax.set_ylim(min_energy, max_energy)

            ani = animation.FuncAnimation(self.fig, self.set_data, frames=timeIdxs, \
                                              fargs=(t,U), interval=self.intervalTime, \
                                              blit=True, repeat=False)

        elif (self.plot_type == 'E_vs_t'):
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Total energy (J)")
            self.ax.set_title("Total Energy vs. Time")

            t = self.pm.t
            K = self.pm.K
            U = self.pm.U
            E = K+U

            min_energy = np.min([np.min(K), np.min(U), np.min(E)])
            max_energy = np.max([np.max(K), np.max(U), np.max(E)])
            min_energy, max_energy = self.findMinMax(min_energy, max_energy)
            
            self.ax.set_xlim(0, np.max(t)*self.SCALE_FACTOR)
            self.ax.set_ylim(min_energy, max_energy)

            ani = animation.FuncAnimation(self.fig, self.set_data, frames=timeIdxs, \
                                              fargs=(t,E), interval=self.intervalTime, \
                                              blit=True, repeat=False)

        return ani


    def findMinMax(self, minVal, maxVal):

        if (minVal < 0 and maxVal < 0):
            minVal = minVal*self.SCALE_FACTOR
            maxVal = maxVal/self.SCALE_FACTOR
        elif (minVal < 0 and maxVal >= 0):
            minVal = minVal*self.SCALE_FACTOR
            maxVal = maxVal*self.SCALE_FACTOR
        elif (minVal >= 0 and maxVal >= 0):
            minVal = minVal/self.SCALE_FACTOR
            maxVal = maxVal*self.SCALE_FACTOR

        return minVal, maxVal
            
        
