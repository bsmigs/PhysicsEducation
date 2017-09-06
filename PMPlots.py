import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class PMPlots:

    def __init__(self, ax, pm):
        self.line, = ax.plot([], [], 'k-')
        self.ax = ax
        self.x = pm.pos[:,0]
        self.y = pm.pos[:,1]
        self.t = pm.t

        self.ax.set_xlim(0, np.max(self.x)*1.05)
        self.ax.set_ylim(0, np.max(self.y)*1.05)
        self.ax.grid(True)
        self.time_text = self.ax.text(0.02, 0.95, '', transform=self.ax.transAxes)


    def set_data(self, i):
        self.line.set_data(self.x[0:i+1], self.y[0:i+1])
        self.time_text.set_text('time = %.3f (s)' % self.t[i])
        return self.line, self.time_text
    

    def set_plot_type(self, plot_type):
        if (plot_type == 'y_vs_x'):
            self.ax.set_xlabel("Horizontal position (m)")
            self.ax.set_ylabel("Vertical position (m)")
            self.ax.set_title("Vertical Position vs. Horizontal Position")
            
        elif (plot_type == 'x_vs_t'):
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Horizontal position (m)")
            self.ax.set_title("Horizontal Position vs. Time")
            
        elif (plot_type == 'y_vs_t'):
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Vertical position (m)")
            self.ax.set_title("Vertical Position vs. Time")

    
        
