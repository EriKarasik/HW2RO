from FandIK import *
from sympy import *

def cross(a, b): return Matrix([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])
def Tzd(): return Matrix([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0]])
def Tyd(): return Matrix([[0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])
def Txd(): return Matrix([[0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
def Rxd(a): return Matrix([[0,0,0,0],[0,-sin(a),-cos(a),0],[0,cos(a),-sin(a),0],[0,0,0,0]])
def Ryd(a): return Matrix([[-sin(a),0,cos(a),0],[0,0,0,0],[-cos(a),0,-sin(a),0],[0,0,0,0]])
def Rzd(a): return Matrix([[-sin(a),-cos(a),0,0],[cos(a),-sin(a),0,0],[0,0,0,0],[0,0,0,0]])
def Rm1(): return Matrix([[cos(q0),sin(q0),0,0],[-sin(q0)*cos(q1),cos(q0)*cos(q1),sin(q1),0], #R^-1
                          [sin(q0)*sin(q1),-cos(q0)*sin(q1),cos(q1),0],[0,0,0,1]])

q0, q1, q2, a1, a2 = symbols('q0 q1 q2 a1 a2')
H=simplify(Rz(q0)*Tz(a1)*Rx(q1)*Ty(a2)*Ty(q2))
Td = Matrix(simplify(Rzd(q0)*Tz(a1)*Rx(q1)*Ty(a2)*Ty(q2)*Rm1()))
J1 = Matrix([Td[0,3], Td[1,3], Td[2,3], Td[2,1], Td[0,2], Td[1,0]])
Td = Matrix(simplify(Rz(q0)*Tz(a1)*Rxd(q1)*Ty(a2)*Ty(q2)*Rm1()))
J2 = Matrix([Td[0,3], Td[1,3], Td[2,3], Td[2,1], Td[0,2], Td[1,0]])
Td = Matrix(simplify(Rz(q0)*Tz(a1)*Rx(q1)*Ty(a2)*Tyd()*Rm1()))
J3 = Matrix([Td[0,3], Td[1,3], Td[2,3], Td[2,1], Td[0,2], Td[1,0]])
Jcb1 = J1.col_insert(1, J2.col_insert(1, J3))

T00 = eye(4)
O0 = Matrix([[T00[0,3]],[T00[1,3]],[T00[2,3]]])
z0 = Matrix([[T00[0,2]],[T00[1,2]],[T00[2,2]]])
T01 = Rz(q0)*Tz(a1)
O1 = Matrix([[T01[0,3]],[T01[1,3]],[T01[2,3]]])
z1 = Matrix([[T01[0,0]],[T01[1,0]],[T01[2,0]]])
T02 = Rz(q0)*Tz(a1)*Rx(q1)*Ty(a2)
O2 = Matrix([[T02[0,3]],[T02[1,3]],[T02[2,3]]])
z2 = Matrix([[T02[0,1]],[T02[1,1]],[T02[2,1]]])
T03 = Rz(q0)*Tz(a1)*Rx(q1)*Ty(a2)*Ty(q2)
O3 = Matrix([[T03[0,3]],[T03[1,3]],[T03[2,3]]])
J1 = simplify(Matrix([cross(z0, O3-O0), z0]))
J2 = simplify(Matrix([cross(z1, O3-O1), z1]))
J3 = simplify(Matrix([z2, 0, 0, 0]))
Jcb2 = J1.col_insert(1, J2.col_insert(1, J3))
# x: (a2 + q2)sq0*cq1
# y: (a2 + q2)cq0*cq1
# z: a1 + (a2 + q2)sq1
Jcb3 = Matrix([
    [-(a2 + q2)*cos(q0)*cos(q1),  (a2 + q2)*sin(q0)*sin(q1), -sin(q0)*cos(q1)],
    [-(a2 + q2)*sin(q0)*cos(q1), -(a2 + q2)*sin(q1)*cos(q0),  cos(q0)*cos(q1)],
    [                         0,          (a2 + q2)*cos(q1),          sin(q1)],
    [                         0,                    cos(q0),                0],
    [                         0,                    sin(q0),                0],
    [                         1,                          0,                0]])
J = Matrix(Transpose(Jcb1)) # rref check only columns, we need to check rows
print(sympify(J.rref())) # to check rows on independency
print(Jcb1-Jcb2, Jcb2-Jcb3) # checking equality of jacobians