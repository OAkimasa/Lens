import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection='3d')

LX = 6
LY = 6
LZ = 6
geneNum = 500
Nair = 1  # 空気の屈折率
Nn = 1.5  # レンズの屈折率
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

# レイトレーシング、光線ベクトルとレンズの交点を持つときの係数Ｔを求める関数
def rayTraceDecideT_Lens0L(startV, directionV):
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

def rayTraceDecideT_Lens0R(startV, directionV):
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
    T = (-B+np.sqrt(B**2-A*C))/A
    return T

# レンズ表面の法線を求める関数
def decideNormalV_Lens0(pointV):
    nornalVx = (2/Rx**2)*pointV[0]
    nornalVy = (2/Ry**2)*pointV[1]
    nornalVz = (2/Rz**2)*pointV[2]
    normalV = np.array([nornalVx, nornalVy, nornalVz])
    return normalV

# スネルの公式から屈折光の方向ベクトルを求める関数
def decideRefractionVL(rayV, normalV, Nin, Nout):
    # 正規化
    rayV = rayV/np.linalg.norm(rayV)
    normalV = normalV/np.linalg.norm(normalV)
    # 係数A
    A = Nin/Nout
    # 入射角
    cos_t_in = -np.dot(rayV,normalV)
    #量子化誤差対策
    if cos_t_in<-1.:
        cos_t_in = -1.
    elif cos_t_in>1.:
        cos_t_in = 1.
    # スネルの法則
    sin_t_in = np.sqrt(1.0 - cos_t_in**2)
    sin_t_out = sin_t_in*A
    if sin_t_out>1.0:
        #全反射する場合
        return np.zeros(3)
    cos_t_out = np.sqrt(1 - sin_t_out**2)
    # 係数B
    B = -cos_t_out + A*cos_t_in
    # 出射光線の方向ベクトル
    outRayV = A*rayV + B*normalV
    # 正規化
    outRayV = outRayV/np.linalg.norm(outRayV)
    return outRayV

def decideRefractionVR(rayV, normalV, Nin, Nout):
    # 正規化
    rayV = rayV/np.linalg.norm(rayV)
    normalV = normalV/np.linalg.norm(normalV)
    # 係数A
    A = Nin/Nout
    # 入射角
    cos_t_in = np.dot(rayV,normalV)
    #量子化誤差対策
    if cos_t_in<-1.:
        cos_t_in = -1.
    elif cos_t_in>1.:
        cos_t_in = 1.
    # スネルの法則
    sin_t_in = np.sqrt(1.0 - cos_t_in**2)
    sin_t_out = sin_t_in*A
    if sin_t_out>1.0:
        #全反射する場合
        return np.zeros(3)
    cos_t_out = np.sqrt(1 - sin_t_out**2)
    # 係数B
    B = -cos_t_out + A*cos_t_in
    # 出射光線の方向ベクトル
    outRayV = A*rayV + B*normalV
    # 正規化
    outRayV = outRayV/np.linalg.norm(outRayV)
    return outRayV

# ２点の位置ベクトルから直線を引く関数
def plotLineRed(startPointV, endPointV):
    startX = startPointV[0]
    startY = startPointV[1]
    startZ = startPointV[2]
    endX = endPointV[0]
    endY = endPointV[1]
    endZ = endPointV[2]
    ax.plot([startX,endX],[startY,endY],[startZ,endZ],
        'o-',ms='2',linewidth=0.5,color='r')

def plotLinePurple(startPointV, endPointV):
    startX = startPointV[0]
    startY = startPointV[1]
    startZ = startPointV[2]
    endX = endPointV[0]
    endY = endPointV[1]
    endZ = endPointV[2]
    ax.plot([startX,endX],[startY,endY],[startZ,endZ],
        'o-',ms='2',linewidth=0.5,color='purple')

def plotLineBlue(startPointV, endPointV):
    startX = startPointV[0]
    startY = startPointV[1]
    startZ = startPointV[2]
    endX = endPointV[0]
    endY = endPointV[1]
    endZ = endPointV[2]
    ax.plot([startX,endX],[startY,endY],[startZ,endZ],
        'o-',ms='2',linewidth=0.5,color='blue')

raySPoint0 = np.array([-6, 2, 0])  # 入射光の始点
directionVector0 = np.array([1, 0, 0])  # 入射光の方向ベクトルを設定
T = rayTraceDecideT_Lens0L(raySPoint0, directionVector0)  # 交点のための係数
rayEPoint0 = raySPoint0 + T*directionVector0  # 入射光の終点
plotLineRed(raySPoint0, rayEPoint0)  # 入射光描画

refractSPoint0 = rayEPoint0  # 入射光の終点を引き継ぐ
normalV0 = decideNormalV_Lens0(refractSPoint0)  # レンズの法線を求める
#ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],normalV0[0],normalV0[1],normalV0[2],linewidth=0.5)
# 屈折光の方向ベクトルを求める
refractionV0 = decideRefractionVL(directionVector0, normalV0, Nair, Nn)
#ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],refractionV0[0],refractionV0[1],refractionV0[2],linewidth=0.5)
# 係数Tを求めて、屈折光の終点も求める
T = rayTraceDecideT_Lens0R(refractSPoint0, refractionV0)
refractEPoint0 = refractSPoint0 + T*refractionV0
plotLinePurple(refractSPoint0,refractEPoint0)  # 屈折光の描画

raySPoint1 = refractEPoint0  # 屈折光の終点を引き継ぐ
normalV1 = decideNormalV_Lens0(raySPoint1)  # レンズの法線を求める
#ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],normalV1[0],normalV1[1],normalV1[2],linewidth=0.5)
# 屈折光の方向ベクトルを求める
refractionV1 = decideRefractionVR(refractionV0, normalV1, Nair, Nn)
#ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],refractionV1[0],refractionV1[1],refractionV1[2],linewidth=0.5)
T = 10  # 衝突判定がないので適当に置いた
rayEPoint1 = raySPoint1 + T*refractionV1  # 屈折光の終点
plotLineBlue(raySPoint1,rayEPoint1)  # 屈折光の描画


# グラフの見た目について
ax.set_xlim(-LX, LX)
ax.set_ylim(-LY, LY)
ax.set_zlim(-LZ, LZ)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.title('Near the Lens')

plt.show()