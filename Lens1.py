import numpy as np
import math
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.widgets import Slider


# 球面レンズを複数重ね、収差の変化を確認する
LX = 4
LY = 4
LZ = 4
geneNum = 500
raysDensity = 0.5  # 入射光線の密度　スライダーにしたい
nn = 1.6  # レンズの屈折率
zplus = 0  # 入射光のz軸方向移動
Rx = 0.5  # 楕円面の倍率
Ry = 3  # 楕円面の倍率
Rz = 3  # 楕円面の倍率

class CalcLine:
    def __init__(self):
        pass

    # 受け取ったx,y,z座標から(x,y,z)の組を作る関数
    def makePoints(self, point0, point1, point2, shape0, shape1):
        result = [None]*(len(point0)+len(point1)+len(point2))
        result[::3] = point0
        result[1::3] = point1
        result[2::3] = point2
        result = np.array(result)
        result = result.reshape(shape0, shape1)
        return result

    # 光線と法線のなす角を計算し、なす角と回転基準ベクトルを返す
    def calcAngle(self,enterRayVector,surfacePointx,surfacePointy,surfacePointz):
        # レンズの法線を計算
        normalVector = []
        normalVectorx = []
        normalVectory = []
        normalVectorz = []
        base = []
        for i in range(len(surfacePointx)):
            normalVectorx.append((2/Rx**2)*surfacePointx[i])
            normalVectory.append((2/Ry**2)*surfacePointy[i])
            normalVectorz.append((2/Rz**2)*surfacePointz[i])
            base.append((2/Rx**2)*surfacePointx[i])
            base.append((2/Ry**2)*surfacePointy[i])
            base.append((2/Rz**2)*surfacePointz[i])
            normalVector.append(base)
            base = []
        normalVectorx = np.array(normalVectorx)
        normalVectory = np.array(normalVectory)
        normalVectorz = np.array(normalVectorz)
        #print(normalVector.shape)

        # 入射光と法線のなす角を計算
        angle = []
        for i in range(len(enterRayVector)):
            angle.append(np.arccos(
                np.dot(enterRayVector[i], normalVector[i])/(
                    np.linalg.norm(enterRayVector[i], ord=2
                    )*np.linalg.norm(normalVector[i], ord=2))))
        # 屈折の法則
        angle = np.arcsin((np.sin(angle)/nn))

        # 回転基準とするベクトルを生成
        Cross = np.cross(enterRayVector, normalVector)
        Norm = []
        for i in Cross:
            Norm.append(np.linalg.norm(i, ord=2))
        rotBaseVector = []
        for i in range(len(Cross)):
            rotBaseVector.append(np.divide(
                Cross[i], Norm[i], out=np.zeros_like(Cross[i]), where=Norm[i] != 0
                ))
        return angle, rotBaseVector, normalVectorx, normalVectory, normalVectorz

    # 入射光の方向ベクトル、回転基準ベクトル、回転角から屈折光を計算する
    def geneRayinLens(self, enterVector, rotBaseVector, rotAngle):
        for i in range(len(enterVector)):
            enterVector[i] = np.array(enterVector[i])
            enterVector[i] = np.append(enterVector[i], 1)
        #print(*enterVector, sep='\n')
        rayinLensBase = []
        rayinLens = []
        for i in range(len(rotAngle)):
            rot44 = [[rotBaseVector[i][0]**2*(1-np.cos(rotAngle[i]))+np.cos(rotAngle[i]),
                    rotBaseVector[i][0]*rotBaseVector[i][1]*(1-np.cos(rotAngle[i]))-rotBaseVector[i][2]*np.sin(rotAngle[i]),
                    rotBaseVector[i][2]*rotBaseVector[i][0]*(1-np.cos(rotAngle[i]))+rotBaseVector[i][1]*np.sin(rotAngle[i]),
                    0],
                    [rotBaseVector[i][0]*rotBaseVector[i][1]*(1-np.cos(rotAngle[i]))+rotBaseVector[i][2]*np.sin(rotAngle[i]),
                    rotBaseVector[i][1]**2*(1-np.cos(rotAngle[i]))+np.cos(rotAngle[i]),
                    rotBaseVector[i][1]*rotBaseVector[i][2]*(1-np.cos(rotAngle[i]))-rotBaseVector[i][0]*np.sin(rotAngle[i]),
                    0],
                    [rotBaseVector[i][2]*rotBaseVector[i][0]*(1-np.cos(rotAngle[i]))-rotBaseVector[i][1]*np.sin(rotAngle[i]),
                    rotBaseVector[i][1]*rotBaseVector[i][2]*(1-np.cos(rotAngle[i]))+rotBaseVector[i][0]*np.sin(rotAngle[i]),
                    rotBaseVector[i][2]**2*(1-np.cos(rotAngle[i]))+np.cos(rotAngle[i]),
                    0],
                    [0,0,0,1]]
            rayinLensBase = np.dot(rot44, enterVector[i])
            rayinLens.append(rayinLensBase)
            rayinLensBase = []
        #print(rayinLens[0])
        rayinLens = np.array(rayinLens)
        #rayinLens = rayinLens*np.linalg.norm(enterVector)  # 厳密性は失われているかも
        #rayinLens = rayinLens*10  # 厳密性は失われているかも
        #print(rayinLens[0])
        return rayinLens


class PlotLens:
    def __init__(self):
        pass

    # 球面レンズが１つの場合
    def plot1Lens(self):
        # 楕円面で球面レンズを再現する
        limitTheta = 2*np.pi  # theta生成数
        limitPhi = np.pi  # phi生成数
        theta = np.linspace(0, limitTheta, geneNum)
        phi = np.linspace(0, limitPhi, geneNum)

        Xs = Rx * np.outer(np.cos(theta), np.sin(phi))
        Ys = Ry * np.outer(np.sin(theta), np.sin(phi))
        Zs = Rz * np.outer(np.ones(np.size(theta)), np.cos(phi))
        ax.plot_wireframe(Xs, Ys, Zs, linewidth=0.2)

        # 光線の生成（入射光
        # 始点を正方形の格子点として、終点をレンズの曲面上とする
        yRayStart, zRayStart = np.meshgrid(
            np.arange(-2, 3, 1), np.arange(-2, 3, 1))
        yRayStart = yRayStart.reshape(25)*raysDensity
        zRayStart = zRayStart.reshape(25)*raysDensity
        xRayStart = np.array([-8]*25)

        yRayEnd, zRayEnd = np.meshgrid(
            np.arange(-2, 3, 1), np.arange(-2, 3, 1))
        yRayEnd = yRayEnd.reshape(25)*raysDensity
        zRayEnd = zRayEnd.reshape(25)*raysDensity + zplus
        xRayEnd = -Rx*np.sqrt(1-(yRayEnd/Ry)**2-(zRayEnd/Rz)**2)

        CL = CalcLine()  # インスタンス化
        raysStart = CL.makePoints(xRayStart,yRayStart,zRayStart,25,3)
        raysEnd = CL.makePoints(xRayEnd,yRayEnd,zRayEnd,25,3)

        # 入射光のプロットと次に使うためのベクトルを生成
        RayVector = []
        base = []
        for i in range(len(raysStart)):
            XX = [raysStart[i,0], raysEnd[i,0]]
            YY = [raysStart[i,1], raysEnd[i,1]]
            ZZ = [raysStart[i,2], raysEnd[i,2]]
            #ax.plot(XX, YY, ZZ, 'o-', color='r', ms='2', linewidth=0.5)
            ax.quiver(raysStart[i,0], raysStart[i,1], raysStart[i,2],
            xRayEnd, 0, 0, linewidth=0.5)
            base.append(raysEnd[i,0]-raysStart[i,0])
            base.append(raysEnd[i,1]-raysStart[i,1])
            base.append(raysEnd[i,2]-raysStart[i,2])
            base = base/(
                (raysEnd[i,0]-raysStart[i,0])**2+(
                    raysEnd[i,1]-raysStart[i,1])**2+(
                        raysEnd[i,2]-raysStart[i,2])**2)
            RayVector.append(base)
            base = []
        #print(RayVector)

        # レンズ中の屈折光をプロット
        # 入射光とレンズのなす角を計算
        reflectAngle = CL.calcAngle(
            RayVector, xRayEnd, yRayEnd, zRayEnd)[0]
        rotBaseVector = CL.calcAngle(
            RayVector, xRayEnd, yRayEnd, zRayEnd)[1]
        surfaceNormalx = CL.calcAngle(
            RayVector, xRayEnd, yRayEnd, zRayEnd)[2]
        surfaceNormaly = CL.calcAngle(
            RayVector, xRayEnd, yRayEnd, zRayEnd)[3]
        surfaceNormalz = CL.calcAngle(
            RayVector, xRayEnd, yRayEnd, zRayEnd)[4]

        # 屈折光の始点
        xinLensStart = xRayEnd  # 入射光の終点を引き継ぐ
        yinLensStart = yRayEnd
        zinLensStart = zRayEnd

        rayinLens = CL.geneRayinLens(RayVector, rotBaseVector, reflectAngle)
        #print(rayinLens[0][3], sep='\n')

        yinLensEnd, zinLensEnd = np.meshgrid(
            np.arange(-2, 3, 1), np.arange(-2, 3, 1))
        yinLensEnd = yinLensEnd.reshape(25)*raysDensity  # 要変更
        zinLensEnd = zinLensEnd.reshape(25)*raysDensity + zplus

        for i in range(len(rayinLens)):
            yinLensEnd[i] = rayinLens[i][1]  # 傾きを追加？
            zinLensEnd[i] = zinLensStart[i] - RayVector[i][2] + rayinLens[i][2]
        xinLensEnd = Rx*np.sqrt(1-(yinLensEnd/Ry)**2-(zinLensEnd/Rz)**2)
        print(zinLensStart)

        raysLensStart = CL.makePoints(
            xinLensStart, yinLensStart, zinLensStart, 25, 3)
        raysLensEnd = CL.makePoints(
            xinLensEnd, yinLensEnd, zinLensEnd, 25, 3)
        # 屈折光のプロット
        RayVector = []
        base = []
        for i in range(len(raysLensStart)):  # 関数にまとめてCalcLineクラスへ
            XX = [raysLensStart[i,0], raysLensEnd[i,0]]
            YY = [raysLensStart[i,1], raysLensEnd[i,1]]
            ZZ = [raysLensStart[i,2], raysLensEnd[i,2]]
            #ax.plot(XX, YY, ZZ, 'o-', color='purple', ms='2', linewidth=0.5)
            ax.quiver(raysLensStart[i,0], raysLensStart[i,1], raysLensStart[i,2],
            rayinLens[i][0], rayinLens[i][1], rayinLens[i][2], linewidth=0.5)
            base.append(raysLensEnd[i,0]-raysLensStart[i,0])
            base.append(raysLensEnd[i,1]-raysLensStart[i,1])
            base.append(raysLensEnd[i,2]-raysLensStart[i,2])
            base = base/(
                (raysLensEnd[i,0]-raysLensStart[i,0])**2+(
                    raysLensEnd[i,1]-raysLensStart[i,1])**2+(
                        raysLensEnd[i,2]-raysLensStart[i,2])**2)
            RayVector.append(base)
            base = []

        # ２段階目の屈折光をプロット
        # １段目屈折光とレンズのなす角を計算
        reflectAngle = CL.calcAngle(
            RayVector, xinLensEnd, yinLensEnd, zinLensEnd)[0]
        rotBaseVector = CL.calcAngle(
            RayVector, xinLensEnd, yinLensEnd, zinLensEnd)[1]
        surfaceNormalx = CL.calcAngle(
            RayVector, xRayEnd, yRayEnd, zRayEnd)[2]
        surfaceNormaly = CL.calcAngle(
            RayVector, xRayEnd, yRayEnd, zRayEnd)[3]
        surfaceNormalz = CL.calcAngle(
            RayVector, xRayEnd, yRayEnd, zRayEnd)[4]

        # ２段階目屈折光の始点
        xoutLensStart = xinLensEnd  # １段階目屈折光の終点を引き継ぐ
        youtLensStart = yinLensEnd
        zoutLensStart = zinLensEnd

        rayoutLens = CL.geneRayinLens(RayVector, rotBaseVector, reflectAngle)

        youtLensEnd, zoutLensEnd = np.meshgrid(
            np.arange(-2, 3, 1), np.arange(-2, 3, 1))
        youtLensEnd = youtLensEnd.reshape(25)*raysDensity  # 要変更
        zoutLensEnd = zoutLensEnd.reshape(25)*raysDensity + zplus

        for i in range(len(rayoutLens)):
            youtLensEnd[i] = -surfaceNormaly[i] + youtLensStart[i] - rayoutLens[i][1]
            zoutLensEnd[i] = -surfaceNormalz[i] + zoutLensStart[i] - rayoutLens[i][2]
        xoutLensEnd = np.array([4]*25)

        raysStart = CL.makePoints(
            xoutLensStart, youtLensStart, zoutLensStart, 25, 3)
        raysEnd = CL.makePoints(
            xoutLensEnd, youtLensEnd, zoutLensEnd, 25, 3)
        # 屈折光のプロット
        RayVector = []
        base = []
        for i in range(len(raysLensStart)):  # 関数にまとめてCalcLineクラスへ
            XX = [raysStart[i,0], raysEnd[i,0]]
            YY = [raysStart[i,1], raysEnd[i,1]]
            ZZ = [raysStart[i,2], raysEnd[i,2]]
            ax.plot(XX, YY, ZZ, 'o-', color='b', ms='2', linewidth=0.5)

        # グラフの見た目について
        ax.set_xlim(-LX, LX)
        ax.set_ylim(-LY, LY)
        ax.set_zlim(-LZ, LZ)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.title('Near the Lens')
        '''
        # 目盛り幅を揃える
        max_range = np.array(
            [X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() * 0.5
        mid_x = (X.max()+X.min()) * 0.5
        mid_y = (Y.max()+Y.min()) * 0.5
        mid_z = (Z.max()+Z.min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        '''


if __name__ == "__main__":
    start = time.time()

    PL = PlotLens()  # インスタンス化
    fig = plt.figure(figsize=(16, 8))

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    PL.plot1Lens()  # レンズ付近に注目したい

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    #PL.()  # 焦点付近に注目したい

    print(time.time()-start)
    # グラフ描画
    plt.show()
