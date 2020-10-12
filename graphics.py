from FandIK import *
from numpy import radians
import matplotlib.pyplot as plt

def Jcb(q0, q1, q2):
    return Matrix([
    [-(a2 + q2)*cos(q0)*cos(q1),  (a2 + q2)*sin(q0)*sin(q1), -sin(q0)*cos(q1)],
    [-(a2 + q2)*sin(q0)*cos(q1), -(a2 + q2)*sin(q1)*cos(q0),  cos(q0)*cos(q1)],
    [                         0,          (a2 + q2)*cos(q1),          sin(q1)],
    [                         0,                    cos(q0),                0],
    [                         0,                    sin(q0),                0],
    [                         1,                          0,                0]])

x, y, z, d, t = [], [], [], [], []
for i in range(3600): #3600 becouse function range can't work with float
    q = Matrix([[cos(radians(i/10))], [-2*sin(radians(2*i/10))], [3*cos(radians(3*i/10))]])
    v = Jcb(sin(radians(i/10)), cos(radians(2 * i/10)), sin(radians(3 * i/10)))*q
    t.append(i/10)
    x.append(v[0])
    y.append(v[1])
    z.append(v[2])
    d.append(sqrt(v[0]**2+v[1]**2+v[2]**2)-20)

fig, axes = plt.subplots(4, 1, figsize=(16, 8));
axes[0].scatter(t, x, s=1)
axes[1].scatter(t, y, s=1)
axes[2].scatter(t, z, s=1)
axes[3].scatter(t, d, s=1)
axes[0].set_ylabel('x speed')
axes[1].set_ylabel('y speed')
axes[2].set_ylabel('z speed')
axes[3].set_ylabel('common speed')
axes[-1].set_xlabel('angle')
plt.show()