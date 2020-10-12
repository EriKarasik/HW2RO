from sympy import Matrix, cos, sin, atan2, sqrt, zeros, pi

class dot:
    def __init__(self, x, y, z): self.x, self.y, self.z = self, y, z
    def getDistance(self, other): return sqrt((self.x-other.x)**2+(self.y-other.y)**2+(self.z-other.z)**2)
    def checkSphere(self): return dot.getDistance(self, dot(0,0,10)) <= a2+dmax

a1, a2, dmax = 10, 20, 15

def Rx(a): return Matrix([[1, 0, 0, 0], [0, cos(a), -sin(a), 0], [0, sin(a), cos(a), 0], [0, 0, 0, 1]])
def Ry(a): return Matrix([[cos(a), 0, sin(a), 0], [0, 1, 0, 0], [-sin(a), 0, cos(a), 0], [0, 0, 0, 1]])
def Rz(a): return Matrix([[cos(a), -sin(a), 0, 0], [sin(a), cos(a), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
def Tz(a): return Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, a], [0, 0, 0, 1]])
def Ty(a): return Matrix([[1, 0, 0, 0], [0, 1, 0, a], [0, 0, 1, 0], [0, 0, 0, 1]])
def Tx(a): return Matrix([[1, 0, 0, a], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
def FK(q): return Rz(q[0])*Tz(a1)*Rx(q[1])*Ty(a2)*Ty(q[2])
def IK(goal):
    q1, q2 = zeros(3), zeros(3)
    if dot.checkSphere(goal):
        q1[0] = atan2(goal.y, goal.x)
        q2[0] = q1[0] + pi # each goal can be reached two ways: elbow up and elbow down
        current = dot(0, 0, a1)
        q1[2] = dot.getDistance(goal, current) - a2
        q2[2] = q1[2]
        q1[1] = atan2(goal.z-a1, sqrt(goal.x**2+goal.y**2))
        q2[1] = pi - q1[1]
    else: print("can't reach the goal")
    return q1, q2