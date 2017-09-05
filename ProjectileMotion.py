"""
General Numerical Solver for 2D-projectile motion.

author: Brian Smigielski
email: bsmigs@gmail.com
website: http://rfground.wordpress.com
"""

import numpy as np
import scipy.integrate as integrate
from scipy import interpolate
import Tkinter as tk

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
        beta = dragCoeff * np.sqrt(dummyState[1]*dummyState[1] + dummyState[3]*dummyState[3])

        # x velocity
        dxdt = dummyState[1]
        # y velocity
        dydt = dummyState[3]
        
        # x-acceleration
        dvxdt = -beta * dummyState[1]
        # y-acceleration
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

        
