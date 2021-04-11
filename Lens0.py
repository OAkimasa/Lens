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
'''
    # 光線を傾きと切片から生成するメソッド
    def generateLine(self):

    # レンズによる屈折角を計算するメソッド
    def calcAngle(self):
'''

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
        zRayEnd = zRayEnd.reshape(25)*raysDensity
        xRayEnd = -Rx*np.sqrt(1-(yRayEnd/Ry)**2-(zRayEnd/Rz)**2)

        CL = CalcLine()  # インスタンス化
        raysStart = CL.makePoints(xRayStart,yRayStart,zRayStart,25,3)
        raysEnd = CL.makePoints(xRayEnd,yRayEnd,zRayEnd,25,3)
        print(raysStart)

        # 入射光のプロット
        for i in range(len(raysStart)):
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
