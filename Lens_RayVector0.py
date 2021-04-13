import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection='3d')

LX = 6
LY = 6
LZ = 6
geneNum = 500
Rx = 0.5  # 楕円面の倍率
Ry = 3  # 楕円面の倍率
Rz = 3  # 楕円面の倍率

# まずは１本の光線で作る
# 楕円面で球面レンズを再現する
limitTheta = 2*np.pi  # theta生成数
limitPhi = np.pi  # phi生成数
theta = np.linspace(0, limitTheta, geneNum)
phi = np.linspace(0, limitPhi, geneNum)

Xs = Rx * np.outer(np.cos(theta), np.sin(phi))
Ys = Ry * np.outer(np.sin(theta), np.sin(phi))
Zs = Rz * np.outer(np.ones(np.size(theta)), np.cos(phi))
ax.plot_wireframe(Xs, Ys, Zs, linewidth=0.2)

raySPoint0 = np.array([-6, 1, 1])  # 入射光の始点
directionVector0 = np.array([1, 0.2, 0])  # 入射光の方向ベクトルを設定

# レイトレーシング、光線ベクトルとレンズの交点を持つときの係数Ｔを求める関数
def rayTraceDecideT(startV, directionV):
    A = (directionV[0]**2/Rx**2)+(
            directionV[1]**2/Ry**2)+(
            directionV[2]**2/Rz**2)
    #print(A)
    B = (startV[0]*directionV[0]/Rx**2)+(
            startV[1]*directionV[1]/Ry**2)+(
            startV[2]*directionV[2]/Rz**2)
    #print(B)
    C = -1+(startV[0]**2/Rx**2)+(
            startV[1]**2/Ry**2)+(
            startV[2]**2/Rz**2)
    #print(C)
    T = (-B-np.sqrt(B**2-A*C))/A
    return T

T = rayTraceDecideT(raySPoint0, directionVector0)
rayEPoint0 = raySPoint0 + T*directionVector0  # 入射光の終点

# ２点の位置ベクトルから直線を引く関数
def plotLine(startPointV, endPointV):
    startX = startPointV[0]
    startY = startPointV[1]
    startZ = startPointV[2]
    endX = endPointV[0]
    endY = endPointV[1]
    endZ = endPointV[2]
    ax.plot([startX,endX],[startY,endY],[startZ,endZ],linewidth=0.5,color='r')

plotLine(raySPoint0, rayEPoint0)  # 入射光描画


# グラフの見た目について
ax.set_xlim(-LX, LX)
ax.set_ylim(-LY, LY)
ax.set_zlim(-LZ, LZ)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.title('Near the Lens')

plt.show()