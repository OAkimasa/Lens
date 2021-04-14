import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import time


LX = 4
LY = 4
LZ = 4
geneNum = 500
Nair = 1  # 空気の屈折率
Nn = 1.5  # レンズの屈折率
Rx = 0.5  # 楕円面の倍率
Ry = 3  # 楕円面の倍率
Rz = 3  # 楕円面の倍率
centerX = -6  # 入射光の始点の中心座標
centerY = 0  # 入射光の始点の中心座標
centerZ = 0  # 入射光の始点の中心座標
rayDensity = 0.3  # 入射光の密度
focusX = 10.75  # 焦点付近の描画範囲を平行移動
lens0V = np.array([0, 0, 0])

class VectorFunctions:
    # 受け取ったx,y,z座標から(x,y,z)の組を作る関数
    def makePoints(self, point0, point1, point2, shape0, shape1):
        result = [None]*(len(point0)+len(point1)+len(point2))
        result[::3] = point0
        result[1::3] = point1
        result[2::3] = point2
        result = np.array(result)
        result = result.reshape(shape0, shape1)
        return result

    # レイトレーシング、光線ベクトルとレンズ0の交点を持つときの係数Ｔを求める関数
    def rayTraceDecideT_Lens0L(self, startV, directionV):
        startV = startV - lens0V
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

    def rayTraceDecideT_Lens0R(self, startV, directionV):
        startV = startV - lens0V
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

    # レンズ0表面の法線を求める関数
    def decideNormalV_Lens0(self, pointV):
        pointV = pointV - lens0V
        nornalVx = (2/Rx**2)*pointV[0]
        nornalVy = (2/Ry**2)*pointV[1]
        nornalVz = (2/Rz**2)*pointV[2]
        normalV = np.array([nornalVx, nornalVy, nornalVz])
        return normalV

    # スネルの公式から屈折光の方向ベクトルを求める関数
    def decideRefractionVL(self, rayV, normalV, Nair, Nn):
        # 正規化
        rayV = rayV/np.linalg.norm(rayV)
        normalV = normalV/np.linalg.norm(normalV)
        # 係数A
        A = Nair/Nn
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

    def decideRefractionVR(self, rayV, normalV, Nair, Nn):
        # 正規化
        rayV = rayV/np.linalg.norm(rayV)
        normalV = normalV/np.linalg.norm(normalV)
        # 係数A
        A = Nair/Nn
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
    def plotLineRed(self, startPointV, endPointV):
        startX = startPointV[0]
        startY = startPointV[1]
        startZ = startPointV[2]
        endX = endPointV[0]
        endY = endPointV[1]
        endZ = endPointV[2]
        ax.plot([startX,endX],[startY,endY],[startZ,endZ],
            'o-',ms='2',linewidth=0.5,color='r')

    def plotLinePurple(self, startPointV, endPointV):
        startX = startPointV[0]
        startY = startPointV[1]
        startZ = startPointV[2]
        endX = endPointV[0]
        endY = endPointV[1]
        endZ = endPointV[2]
        ax.plot([startX,endX],[startY,endY],[startZ,endZ],
            'o-',ms='2',linewidth=0.5,color='purple')

    def plotLineBlue(self, startPointV, endPointV):
        startX = startPointV[0]
        startY = startPointV[1]
        startZ = startPointV[2]
        endX = endPointV[0]
        endY = endPointV[1]
        endZ = endPointV[2]
        ax.plot([startX,endX],[startY,endY],[startZ,endZ],
            'o-',ms='2',linewidth=0.5,color='blue')


class GenerateRays:
    # 第０レンズと光線の描画
    def generateRayLens0(self):
        # 楕円面で球面レンズを再現する
        limitTheta = 2*np.pi  # theta生成数
        limitPhi = np.pi  # phi生成数
        theta = np.linspace(0, limitTheta, geneNum)
        phi = np.linspace(0, limitPhi, geneNum)
        Xs = Rx * np.outer(np.cos(theta), np.sin(phi)) + lens0V[0]
        Ys = Ry * np.outer(np.sin(theta), np.sin(phi)) + lens0V[1]
        Zs = Rz * np.outer(np.ones(np.size(theta)), np.cos(phi)) + lens0V[2]
        ax.plot_wireframe(Xs, Ys, Zs, linewidth=0.2)

        VF = VectorFunctions()  # インスタンス化

        # 始点を生成する
        pointsX = np.array([centerX]*25) + lens0V[0]
        pointsY, pointsZ = np.meshgrid(
            np.arange(-2+centerY, 3+centerY, 1),
            np.arange(-2+centerZ, 3+centerZ, 1))
        pointsY = pointsY.reshape(25)*rayDensity + lens0V[1]
        pointsZ = pointsZ.reshape(25)*rayDensity + lens0V[2]
        raySPoint0 = VF.makePoints(pointsX, pointsY, pointsZ, 25, 3)

        for i in raySPoint0:
            raySPoint0 = i
            directionVector0 = np.array([1, 0, 0])  # 入射光の方向ベクトルを設定
            T = VF.rayTraceDecideT_Lens0L(raySPoint0, directionVector0)  # 交点のための係数
            rayEPoint0 = raySPoint0 + T*directionVector0  # 入射光の終点
            VF.plotLineRed(raySPoint0, rayEPoint0)  # 入射光描画

            refractSPoint0 = rayEPoint0  # 入射光の終点を引き継ぐ
            normalV0 = VF.decideNormalV_Lens0(refractSPoint0)  # レンズの法線を求める
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],normalV0[0],normalV0[1],normalV0[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV0 = VF.decideRefractionVL(directionVector0, normalV0, Nair, Nn)
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],refractionV0[0],refractionV0[1],refractionV0[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens0R(refractSPoint0, refractionV0)
            refractEPoint0 = refractSPoint0 + T*refractionV0
            VF.plotLinePurple(refractSPoint0,refractEPoint0)  # 屈折光の描画

            raySPoint1 = refractEPoint0  # 屈折光の終点を引き継ぐ
            normalV1 = VF.decideNormalV_Lens0(raySPoint1)  # レンズの法線を求める
            #ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],normalV1[0],normalV1[1],normalV1[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV1 = VF.decideRefractionVR(refractionV0, normalV1, Nair, Nn)
            #ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],refractionV1[0],refractionV1[1],refractionV1[2],linewidth=0.5)
            T = 12  # 衝突判定がないので適当に置いた
            rayEPoint1 = raySPoint1 + T*refractionV1  # 屈折光の終点
            VF.plotLineBlue(raySPoint1,rayEPoint1)  # 屈折光の描画

        # グラフの見た目について
        ax.set_xlim(-LX, LX)
        ax.set_ylim(-LY, LY)
        ax.set_zlim(-LZ, LZ)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.title('Near the Lens')

    # 焦点付近を描画
    def focusFocus(self, func):
        func()
        # グラフの見た目について
        ax.set_xlim(-LX/4+lens0V[0]+focusX, LX/4+lens0V[0]+focusX)
        ax.set_ylim(-LY/4+lens0V[1], LY/4+lens0V[1])
        ax.set_zlim(-LZ/4+lens0V[2], LZ/4+lens0V[2])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.title('Focus')


if __name__ == "__main__":
    start = time.time()

    GR = GenerateRays()  # インスタンス化
    fig = plt.figure(figsize=(16, 8))

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    GR.generateRayLens0()  # レンズ付近を描画

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    GR.focusFocus(GR.generateRayLens0)  # 焦点付近を描画

    print(time.time()-start)
    # グラフ描画
    plt.show()