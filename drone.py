"""
@author: Diego Ferruzzo Correa
@date: 07/07/2022
"""
from numpy import array
#from sympy import *


class drone:
    def __init__(self):

        from numpy import zeros, array
        self.Ir = 1e-3      # motor's moment of inertia
        self.Ixx = 16.83e-3  # x-axis inertia
        self.Iyy = 16.38e-3  # y-axis inertia
        self.Izz = 28.34e-3  # z-axis inertia
        self.Omr = 1.0       # propeller's relative angular velocity
        self.l = 1.0       # distance from the rotor to the CG
        self.g = 9.8       # gravity constant
        # simulation parameters
        self.t0 = 0.0  # simulation initial time
        self.tf = 10.0  # simulation end time
        self.dt = 1e-3  # simulation step
        self.x0 = zeros((12, 1))    # initial condition vector

    def dWdt(self, eta, deta):
        # returns dW/dt
        from numpy import array, sin, cos, tan
        phi = eta[0]
        theta = eta[1]
        psi = eta[2]
        dphi = deta[0]
        dtheta = deta[1]
        dpsi = deta[2]
        return array([[0, ((sin(phi)*dtheta/cos(theta)**2)+cos(phi)*tan(theta)*dphi).item(), (-sin(phi)*tan(theta)*dphi+(cos(phi)*dtheta/cos(theta)**2)).item()],
                      [0, (-sin(phi)*dphi).item(), (-cos(phi)*dphi).item()],
                      [0, ((sin(phi)*tan(theta)*dtheta + cos(phi)*dphi)/cos(theta)).item(), ((-sin(phi)*dphi+cos(phi)*tan(theta)*dtheta)/cos(theta)).item()]]).reshape(3, 3)

    def W(self, eta):
        # returns W matrix
        from numpy import array, tan, sin, cos
        phi = eta[0]
        theta = eta[1]
        psi = eta[2]
        return array([[1, (tan(theta)*sin(phi)).item(), (tan(theta)*cos(phi)).item()],
                      [0, cos(phi).item(), -sin(phi).item()],
                      [0, (sin(phi)/cos(theta)).item(), (cos(phi)/cos(theta)).item()]]).reshape(3, 3)

    def dynamics(self, _x):
        from numpy import sin, cos, array, matmul, tan
        from numpy.linalg import inv
        # state vector
        x = _x[0]
        dx = _x[1]
        y = _x[2]
        dy = _x[3]
        z = _x[4]
        dz = _x[5]
        phi = _x[6]
        dphi = _x[7]
        theta = _x[8]
        dtheta = _x[9]
        psi = _x[10]
        dpsi = _x[11]
        #
        xi = array([x, y, z])
        dxi = array([dx, dy, dz])
        eta = array([phi, theta, psi])
        deta = array([dphi, dtheta, dpsi])
        #
        Ir = self.Ir
        Ixx = self.Ixx
        Iyy = self.Iyy
        Izz = self.Izz
        Omr = self.Omr
        #
        a1 = (tan(theta)/(Iyy*Izz))*((Ixx*Iyy-Iyy**2-Ixx*Izz+Izz**2)
                                     * cos(phi)**2+(Ixx+Iyy)*Izz-Izz**2).item()
        a2 = (1/(Iyy*Izz))*(Ixx*Iyy-Iyy**2-Ixx*Izz+Izz**2) * \
            cos(phi)*sin(phi)*sin(theta).item()
        a3 = ((1/cos(theta))+tan(theta)*sin(theta)*sin(phi)**2*((Izz-Ixx)/Iyy)-tan(theta)*sin(theta)
              * cos(phi)**2*((Ixx-Iyy)/Izz)+cos(theta)*(sin(phi)**2-cos(phi)**2)*((Izz-Iyy)/Ixx)).item()
        a4 = (-(Iyy-Izz)*cos(phi)*sin(phi)/(Ixx)).item()
        a5 = (((Iyy**2*Izz-Iyy*Izz**2)*cos(phi)*cos(theta)**2*sin(phi))/(Ixx*Iyy*Izz)-((Ixx**2*Iyy -
              Ixx*Iyy**2-Ixx**2*Izz+Ixx*Izz**2)*cos(phi)*sin(phi)*sin(theta)**2)/(Ixx*Iyy*Izz)).item()
        a6 = (-((Ixx*Iyy-Iyy**2-Ixx*Izz+Izz**2)*cos(phi)*sin(phi))/Iyy*Izz).item()
        a7 = (sin(phi)**2*((Iyy-Ixx)/Izz)+cos(phi) **
              2*((Izz-Ixx)/Iyy)-1)*cos(theta).item()
        a8 = ((Ixx*Iyy-Iyy**2-Ixx*Izz+Izz**2)*cos(phi)
              * sin(phi)*sin(theta)/(Iyy*Izz)).item()
        a9 = ((-(Iyy-Ixx)/Izz)*sin(phi)**2*cos(theta)*sin(theta) -
              ((Izz-Ixx)/Iyy)*cos(phi)**2*cos(theta)*sin(theta)).item()
        a10 = (((Ixx*Iyy - Iyy**2 - Ixx*Izz + Izz**2)*cos(phi)**2 +
               (Ixx + Iyy)*Izz - Izz**2)/(Iyy*Izz*cos(theta))).item()
        a11 = (Ixx*Iyy - Iyy**2 - Ixx*Izz + Izz**2)*cos(phi)*sin(phi)/(Iyy*Izz)
        a12 = -((Ixx*Iyy - Iyy**2 - Ixx*Izz + Izz**2)*cos(phi)**2 +
                (Ixx - Iyy)*Izz - Izz**2)*sin(theta)/(Iyy*Izz*cos(theta))

      # Inertia Matrix
        pass


"""
x0 = array([[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])
mydrone = drone()
print("x0=", x0)
print(mydrone.dynamics(x0))

#print(getattr(mydrone, 'Ir'))
#setattr(mydrone, 'Ir', 2e-3)
#print(getattr(mydrone, 'Ir'))
#print(mydrone.showX0())
#print(zeros((12,1)))
#print(getattr(mydrone, 'x0'))
#print(mydrone.dynamics(zeros((6,1))))

# Reproduzindo o modelo a partir das equações, utilizamos sympy
# os ângulos de Euler

t = Symbol('t')
phi = Function('phi')(t)
theta = Function('theta')(t)
psi = Function('psi')(t)
W = Matrix([[1, tan(theta)*sin(phi), tan(theta)*cos(phi)],[0, cos(phi), -sin(phi)],[0, sin(phi)/cos(theta), cos(phi)/cos(theta)]])
dWdt1 = simplify(Derivative(W, t).doit())
print("W=", pretty(W))
print("dW/dt=\n")
print(pretty(dWdt1))
print(dWdt(1,1,1,1,1,1))
"""
