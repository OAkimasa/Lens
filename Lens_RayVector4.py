# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import time


LX = 4
LY = 4
LZ = 4
geneNum = 500
Nair = 1  # 空気の屈折率
Nn = 1.4  # レンズの屈折率
Nn2 = 1.8
centerX = -8  # 入射光の始点の中心座標
centerY = 0  # 入射光の始点の中心座標
centerZ = 0  # 入射光の始点の中心座標
rayDensity = 0.5  # 入射光の密度
focusX = 4.75  # 焦点付近の描画範囲を平行移動
Rx0 = 0.5  # レンズ０の倍率
Ry0 = 3  # レンズ０の倍率
Rz0 = 3  # レンズ０の倍率
lens0V = np.array([0, 0, 0])  # レンズ０の位置ベクトル
Rx1 = 1  # レンズ１の倍率
Ry1 = 3  # レンズ１の倍率
Rz1 = 3  # レンズ１の倍率
lens1V = np.array([-1.8, 0, 0])  # レンズ１の位置ベクトル
Rx21 = 0.2  # レンズ２の倍率１
Rx22 = 0.2  # レンズ２の倍率２
Ry21 = 3  # レンズ２の倍率１
Ry22 = 2.5  # レンズ２の倍率２
Rz21 = 3  # レンズ２の倍率１
Rz22 = 2.5  # レンズ２の倍率２
lens2V = np.array([-1.1, 0, 0])  # レンズ２の位置ベクトル
Rx31 = 0.1  # レンズ３の倍率１
Rx32 = 0.7  # レンズ３の倍率２
Ry31 = 2.5  # レンズ３の倍率１
Ry32 = 3  # レンズ３の倍率２
Rz31 = 2.5  # レンズ３の倍率１
Rz32 = 3  # レンズ３の倍率２
lens3V = np.array([1.5, 0, 0])  # レンズ３の位置ベクトル
Rx4 = 0.7  # レンズ４の倍率
Ry4 = 3  # レンズ４の倍率
Rz4 = 3  # レンズ４の倍率
lens4V = np.array([1.5, 0, 0])  # レンズ４の位置ベクトル

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
        A = (directionV[0]**2/Rx0**2)+(
                directionV[1]**2/Ry0**2)+(
                directionV[2]**2/Rz0**2)
        #print(A)
        B = (startV[0]*directionV[0]/Rx0**2)+(
                startV[1]*directionV[1]/Ry0**2)+(
                startV[2]*directionV[2]/Rz0**2)
        #print(B)
        C = -1+(startV[0]**2/Rx0**2)+(
                startV[1]**2/Ry0**2)+(
                startV[2]**2/Rz0**2)
        #print(C)
        T = (-B-np.sqrt(B**2-A*C))/A
        return T

    def rayTraceDecideT_Lens0R(self, startV, directionV):
        startV = startV - lens0V
        A = (directionV[0]**2/Rx0**2)+(
                directionV[1]**2/Ry0**2)+(
                directionV[2]**2/Rz0**2)
        #print(A)
        B = (startV[0]*directionV[0]/Rx0**2)+(
                startV[1]*directionV[1]/Ry0**2)+(
                startV[2]*directionV[2]/Rz0**2)
        #print(B)
        C = -1+(startV[0]**2/Rx0**2)+(
                startV[1]**2/Ry0**2)+(
                startV[2]**2/Rz0**2)
        #print(C)
        T = (-B+np.sqrt(B**2-A*C))/A
        return T

    # レンズ0表面の法線を求める関数
    def decideNormalV_Lens0(self, pointV):
        pointV = pointV - lens0V
        nornalVx = (2/Rx0**2)*pointV[0]
        nornalVy = (2/Ry0**2)*pointV[1]
        nornalVz = (2/Rz0**2)*pointV[2]
        normalV = np.array([nornalVx, nornalVy, nornalVz])
        return normalV


    # レイトレーシング、光線ベクトルとレンズ１の交点を持つときの係数Ｔを求める関数
    def rayTraceDecideT_Lens1L(self, startV, directionV):
        startV = startV - lens1V
        A = (directionV[0]**2/Rx1**2)+(
                directionV[1]**2/Ry1**2)+(
                directionV[2]**2/Rz1**2)
        #print(A)
        B = (startV[0]*directionV[0]/Rx1**2)+(
                startV[1]*directionV[1]/Ry1**2)+(
                startV[2]*directionV[2]/Rz1**2)
        #print(B)
        C = -1+(startV[0]**2/Rx1**2)+(
                startV[1]**2/Ry1**2)+(
                startV[2]**2/Rz1**2)
        #print(C)
        T = (-B-np.sqrt(B**2-A*C))/A
        return T

    def rayTraceDecideT_Lens1R(self, startV, directionV):
        T = (lens1V[0] - startV[0])/directionV[0]
        return T

    # レンズ１表面の法線を求める関数
    def decideNormalV_Lens1L(self, pointV):
        pointV = pointV - lens1V
        nornalVx = (2/Rx1**2)*pointV[0]
        nornalVy = (2/Ry1**2)*pointV[1]
        nornalVz = (2/Rz1**2)*pointV[2]
        normalV = np.array([nornalVx, nornalVy, nornalVz])
        return normalV

    def decideNormalV_Lens1R(self, pointV):
        pointV = pointV - lens1V
        nornalVx = 1
        nornalVy = 0
        nornalVz = 0
        normalV = np.array([nornalVx, nornalVy, nornalVz])
        return normalV


    # レイトレーシング、光線ベクトルとレンズ２の交点を持つときの係数Ｔを求める関数
    def rayTraceDecideT_Lens2L(self, startV, directionV):
        startV = startV - lens2V
        A = (directionV[0]**2/Rx21**2)+(
                directionV[1]**2/Ry21**2)+(
                directionV[2]**2/Rz21**2)
        #print(A)
        B = (startV[0]*directionV[0]/Rx21**2)+(
                startV[1]*directionV[1]/Ry21**2)+(
                startV[2]*directionV[2]/Rz21**2)
        #print(B)
        C = -1+(startV[0]**2/Rx21**2)+(
                startV[1]**2/Ry21**2)+(
                startV[2]**2/Rz21**2)
        #print(C)
        T = (-B+np.sqrt(B**2-A*C))/A
        return T

    def rayTraceDecideT_Lens2R(self, startV, directionV):
        startV = startV - lens2V
        A = (directionV[0]**2/Rx22**2)+(
                directionV[1]**2/Ry22**2)+(
                directionV[2]**2/Rz22**2)
        #print(A)
        B = ((startV[0] - 0.8)*directionV[0]/Rx22**2)+(
                startV[1]*directionV[1]/Ry22**2)+(
                startV[2]*directionV[2]/Rz22**2)
        #print(B)
        C = -1+((startV[0] - 0.8)**2/Rx22**2)+(
                startV[1]**2/Ry22**2)+(
                startV[2]**2/Rz22**2)
        #print(C)
        T = (-B-np.sqrt(B**2-A*C))/A
        return T

    # レンズ２表面の法線を求める関数
    def decideNormalV_Lens2L(self, pointV):
        pointV = pointV - lens2V
        nornalVx = -(2/Rx21**2)*pointV[0]
        nornalVy = -(2/Ry21**2)*pointV[1]
        nornalVz = -(2/Rz21**2)*pointV[2]
        normalV = np.array([nornalVx, nornalVy, nornalVz])
        return normalV

    def decideNormalV_Lens2R(self, pointV):
        pointV = pointV - lens2V
        nornalVx = -(2/Rx22**2)*(pointV[0] - 0.8)
        nornalVy = -(2/Ry22**2)*pointV[1]
        nornalVz = -(2/Rz22**2)*pointV[2]
        normalV = np.array([nornalVx, nornalVy, nornalVz])
        return normalV


    # レイトレーシング、光線ベクトルとレンズ３の交点を持つときの係数Ｔを求める関数
    def rayTraceDecideT_Lens3L(self, startV, directionV):
        startV = startV - lens3V
        A = (directionV[0]**2/Rx31**2)+(
                directionV[1]**2/Ry31**2)+(
                directionV[2]**2/Rz31**2)
        #print(A)
        B = ((startV[0] + 1.1)*directionV[0]/Rx31**2)+(
                startV[1]*directionV[1]/Ry31**2)+(
                startV[2]*directionV[2]/Rz31**2)
        #print(B)
        C = -1+((startV[0] + 1.1)**2/Rx31**2)+(
                startV[1]**2/Ry31**2)+(
                startV[2]**2/Rz31**2)
        #print(C)
        T = (-B+np.sqrt(B**2-A*C))/A
        return T

    def rayTraceDecideT_Lens3R(self, startV, directionV):
        startV = startV - lens3V
        A = (directionV[0]**2/Rx32**2)+(
                directionV[1]**2/Ry32**2)+(
                directionV[2]**2/Rz32**2)
        #print(A)
        B = (startV[0]*directionV[0]/Rx32**2)+(
                startV[1]*directionV[1]/Ry32**2)+(
                startV[2]*directionV[2]/Rz32**2)
        #print(B)
        C = -1+(startV[0]**2/Rx32**2)+(
                startV[1]**2/Ry32**2)+(
                startV[2]**2/Rz32**2)
        #print(C)
        T = (-B-np.sqrt(B**2-A*C))/A
        return T

    # レンズ２表面の法線を求める関数
    def decideNormalV_Lens3L(self, pointV):
        pointV = pointV - lens3V
        nornalVx = -(2/Rx31**2)*(pointV[0] + 1.6)
        nornalVy = -(2/Ry31**2)*pointV[1]
        nornalVz = -(2/Rz31**2)*pointV[2]
        normalV = np.array([nornalVx, nornalVy, nornalVz])
        return normalV

    def decideNormalV_Lens3R(self, pointV):
        pointV = pointV - lens3V
        nornalVx = -(2/Rx32**2)*pointV[0]
        nornalVy = -(2/Ry32**2)*pointV[1]
        nornalVz = -(2/Rz32**2)*pointV[2]
        normalV = np.array([nornalVx, nornalVy, nornalVz])
        return normalV


    # レイトレーシング、光線ベクトルとレンズ４の交点を持つときの係数Ｔを求める関数
    def rayTraceDecideT_Lens4L(self, startV, directionV):
        startV = startV - lens4V
        A = (directionV[0]**2/Rx4**2)+(
                directionV[1]**2/Ry4**2)+(
                directionV[2]**2/Rz4**2)
        #print(A)
        B = (startV[0]*directionV[0]/Rx4**2)+(
                startV[1]*directionV[1]/Ry4**2)+(
                startV[2]*directionV[2]/Rz4**2)
        #print(B)
        C = -1+(startV[0]**2/Rx4**2)+(
                startV[1]**2/Ry4**2)+(
                startV[2]**2/Rz4**2)
        #print(C)
        T = (-B-np.sqrt(B**2-A*C))/A
        return T

    def rayTraceDecideT_Lens4R(self, startV, directionV):
        startV = startV - lens4V
        A = (directionV[0]**2/Rx4**2)+(
                directionV[1]**2/Ry4**2)+(
                directionV[2]**2/Rz4**2)
        #print(A)
        B = (startV[0]*directionV[0]/Rx4**2)+(
                startV[1]*directionV[1]/Ry4**2)+(
                startV[2]*directionV[2]/Rz4**2)
        #print(B)
        C = -1+(startV[0]**2/Rx4**2)+(
                startV[1]**2/Ry4**2)+(
                startV[2]**2/Rz4**2)
        #print(C)
        T = (-B+np.sqrt(B**2-A*C))/A
        return T

    # レンズ４表面の法線を求める関数
    def decideNormalV_Lens4L(self, pointV):
        pointV = pointV - lens4V
        nornalVx = (2/Rx4**2)*pointV[0]
        nornalVy = (2/Ry4**2)*pointV[1]
        nornalVz = (2/Rz4**2)*pointV[2]
        normalV = np.array([nornalVx, nornalVy, nornalVz])
        return normalV

    def decideNormalV_Lens4R(self, pointV):
        pointV = pointV - lens4V
        nornalVx = (2/Rx4**2)*pointV[0]
        nornalVy = (2/Ry4**2)*pointV[1]
        nornalVz = (2/Rz4**2)*pointV[2]
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
        Xs = Rx0 * np.outer(np.cos(theta), np.sin(phi)) + lens0V[0]
        Ys = Ry0 * np.outer(np.sin(theta), np.sin(phi)) + lens0V[1]
        Zs = Rz0 * np.outer(np.ones(np.size(theta)), np.cos(phi)) + lens0V[2]
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

    # トリプレット（テッサー）と光線の描画
    def generateRayTriplet(self):
        def plotLensTriplet():
            # １枚目の凸レンズを再現する
            limitTheta = 2*np.pi  # theta生成数
            limitPhi = np.pi  # phi生成数
            theta = np.linspace(0, limitTheta, geneNum)
            phi = np.linspace(0, limitPhi, geneNum)
            Xs = Rx1 * np.outer(np.cos(theta), np.sin(phi))
            Xs = np.where(0<Xs, 0, Xs) + lens1V[0]
            Ys = Ry1 * np.outer(np.sin(theta), np.sin(phi)) + lens1V[1]
            Zs = Rz1 * np.outer(np.ones(np.size(theta)), np.cos(phi)) + lens1V[2]
            ax.plot_wireframe(Xs, Ys, Zs, linewidth=0.2)

            # ２枚目の凹レンズを再現する
            limitTheta = 2*np.pi  # theta生成数
            limitPhi = np.pi  # phi生成数
            theta = np.linspace(0, limitTheta, geneNum)
            phi = np.linspace(0, limitPhi, geneNum)
            Xs = Rx21 * np.outer(np.cos(theta), np.sin(phi))
            Ys = Ry21 * np.outer(np.sin(theta), np.sin(phi))
            Zs = Rz21 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Xs1 = Rx21 * np.outer(np.cos(theta), np.sin(phi))
            Xs2 = Rx22 * np.outer(np.cos(theta), np.sin(phi)) + 0.8
            Ys1 = Ry21 * np.outer(np.sin(theta), np.sin(phi))
            Ys2 = Ry22 * np.outer(np.sin(theta), np.sin(phi))
            Zs1 = Rz21 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Zs2 = Rz22 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Ys = np.where(Xs<0, Ys2, Ys1) + lens2V[1]
            Zs = np.where(Xs<0, Zs2, Zs1) + lens2V[2]
            Xs = np.where(Xs<0, Xs2, Xs1) + lens2V[0]
            ax.plot_wireframe(Xs, Ys, Zs, linewidth=0.2)

            # ３枚目の凹レンズを再現する
            limitTheta = 2*np.pi  # theta生成数
            limitPhi = np.pi  # phi生成数
            theta = np.linspace(0, limitTheta, geneNum)
            phi = np.linspace(0, limitPhi, geneNum)
            Xs = Rx31 * np.outer(np.cos(theta), np.sin(phi))
            Ys = Ry31 * np.outer(np.sin(theta), np.sin(phi))
            Zs = Rz31 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Xs1 = Rx31 * np.outer(np.cos(theta), np.sin(phi))
            Xs2 = Rx32 * np.outer(np.cos(theta), np.sin(phi)) + 1.1
            Ys1 = Ry31 * np.outer(np.sin(theta), np.sin(phi))
            Ys2 = Ry32 * np.outer(np.sin(theta), np.sin(phi))
            Zs1 = Rz31 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Zs2 = Rz32 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Ys = np.where(Xs<0, Ys2, Ys1) + lens3V[1]
            Zs = np.where(Xs<0, Zs2, Zs1) + lens3V[2]
            Xs = np.where(Xs<0, Xs2, Xs1) + lens3V[0] - 1.1
            ax.plot_wireframe(Xs, Ys, Zs, linewidth=0.2)

            # ４枚目のレンズを再現する
            limitTheta = 2*np.pi  # theta生成数
            limitPhi = np.pi  # phi生成数
            theta = np.linspace(0, limitTheta, geneNum)
            phi = np.linspace(0, limitPhi, geneNum)
            Xs = Rx4 * np.outer(np.cos(theta), np.sin(phi)) + lens4V[0]
            Ys = Ry4 * np.outer(np.sin(theta), np.sin(phi)) + lens4V[1]
            Zs = Rz4 * np.outer(np.ones(np.size(theta)), np.cos(phi)) + lens4V[2]
            ax.plot_wireframe(Xs, Ys, Zs, linewidth=0.2)

        plotLensTriplet()
        VF = VectorFunctions()  # インスタンス化

        # 始点を生成する
        pointsX = np.array([centerX]*81) + lens1V[0]
        pointsY, pointsZ = np.meshgrid(
            np.arange(-4+centerY, 5+centerY, 1),
            np.arange(-4+centerZ, 5+centerZ, 1))
        pointsY = pointsY.reshape(81)*rayDensity + lens1V[1]
        pointsZ = pointsZ.reshape(81)*rayDensity + lens1V[2]
        raySPoint0 = VF.makePoints(pointsX, pointsY, pointsZ, 81, 3)

        for i in raySPoint0:
            raySPoint0 = i
            directionVector0 = np.array([1, 0, 0])  # 入射光の方向ベクトルを設定
            T = VF.rayTraceDecideT_Lens1L(raySPoint0, directionVector0)  # 交点のための係数
            rayEPoint0 = raySPoint0 + T*directionVector0  # 入射光の終点
            VF.plotLineRed(raySPoint0, rayEPoint0)  # 入射光描画
            refractSPoint0 = rayEPoint0  # 入射光の終点を引き継ぐ。以下レンズ１についての計算
            normalV_Lens1L = VF.decideNormalV_Lens1L(refractSPoint0)  # レンズの法線を求める
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],normalV_Lens1L[0],normalV_Lens1L[1],normalV_Lens1L[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens1L = VF.decideRefractionVL(directionVector0, normalV_Lens1L, Nair, Nn)
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],refractionV_Lens1L[0],refractionV_Lens1L[1],refractionV_Lens1L[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens1R(refractSPoint0, refractionV_Lens1L)
            refractEPoint0 = refractSPoint0 + T*refractionV_Lens1L
            VF.plotLinePurple(refractSPoint0,refractEPoint0)  # 屈折光の描画
            raySPoint1 = refractEPoint0  # 屈折光の終点を引き継ぐ
            normalV1 = VF.decideNormalV_Lens1R(raySPoint1)  # レンズの法線を求める
            #ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],normalV1[0],normalV1[1],normalV1[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens1R = VF.decideRefractionVR(refractionV_Lens1L, normalV1, Nair, Nn)
            #ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],refractionV_Lens1R[0],refractionV_Lens1R[1],refractionV_Lens1R[2],linewidth=0.5)
            T = VF.rayTraceDecideT_Lens2L(raySPoint1, refractionV_Lens1R)
            rayEPoint1 = raySPoint1 + T*refractionV_Lens1R  # 空気中の屈折光の終点
            VF.plotLineBlue(raySPoint1,rayEPoint1)  # 空気中の屈折光の描画

            refractSPoint_Lens2L = rayEPoint1  # 以下、レンズ２についての計算
            normalV_Lens2L = VF.decideNormalV_Lens2L(refractSPoint_Lens2L)  # レンズの法線を求める
            #ax.quiver(rayEPoint1[0],rayEPoint1[1],rayEPoint1[2],normalV_Lens2L[0],normalV_Lens2L[1],normalV_Lens2L[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens2L = VF.decideRefractionVL(refractionV_Lens1R, normalV_Lens2L, Nair, Nn)
            #ax.quiver(rayEPoint1[0],rayEPoint1[1],rayEPoint1[2],refractionV_Lens2L[0],refractionV_Lens2L[1],refractionV_Lens2L[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens2R(refractSPoint_Lens2L, refractionV_Lens2L)
            refractEPoint_Lens2L = refractSPoint_Lens2L + T*refractionV_Lens2L
            VF.plotLinePurple(refractSPoint_Lens2L,refractEPoint_Lens2L)  # 屈折光の描画
            raySPoint_Lens2R = refractEPoint_Lens2L
            normalV_Lens2R = VF.decideNormalV_Lens2R(raySPoint_Lens2R)
            #ax.quiver(raySPoint_Lens2R[0],raySPoint_Lens2R[1],raySPoint_Lens2R[2],normalV_Lens2R[0],normalV_Lens2R[1],normalV_Lens2R[2],linewidth=0.5)
            refractionV_Lens2R = VF.decideRefractionVR(refractionV_Lens2L, normalV_Lens2R, Nair, Nn)
            #ax.quiver(raySPoint_Lens2R[0],raySPoint_Lens2R[1],raySPoint_Lens2R[2],refractionV_Lens2R[0],refractionV_Lens2R[1],refractionV_Lens2R[2],linewidth=0.5)
            T = VF.rayTraceDecideT_Lens3L(raySPoint_Lens2R, refractionV_Lens2R)
            rayEPoint_Lens3L = raySPoint_Lens2R + T*refractionV_Lens2R
            VF.plotLineRed(raySPoint_Lens2R, rayEPoint_Lens3L)

            refractSPoint_Lens3L = rayEPoint_Lens3L  # 以下、レンズ３についての計算
            normalV_Lens3L = VF.decideNormalV_Lens3L(refractSPoint_Lens3L)  # レンズの法線を求める
            #ax.quiver(refractSPoint_Lens3L[0],refractSPoint_Lens3L[1],refractSPoint_Lens3L[2],normalV_Lens3L[0],normalV_Lens3L[1],normalV_Lens3L[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens3L = VF.decideRefractionVL(refractionV_Lens2R, normalV_Lens3L, Nair, Nn)
            #ax.quiver(refractSPoint_Lens3L[0],refractSPoint_Lens3L[1],refractSPoint_Lens3L[2],refractionV_Lens3L[0],refractionV_Lens3L[1],refractionV_Lens3L[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens3R(refractSPoint_Lens3L, refractionV_Lens3L)
            refractEPoint_Lens3L = refractSPoint_Lens3L + T*refractionV_Lens3L
            VF.plotLinePurple(refractSPoint_Lens3L,refractEPoint_Lens3L)  # 屈折光の描画
            raySPoint_Lens3R = refractEPoint_Lens3L
            normalV_Lens3R = VF.decideNormalV_Lens3R(raySPoint_Lens3R)
            #ax.quiver(raySPoint_Lens3R[0],raySPoint_Lens3R[1],raySPoint_Lens3R[2],normalV_Lens3R[0],normalV_Lens3R[1],normalV_Lens3R[2],linewidth=0.5)
            refractionV_Lens3R = VF.decideRefractionVR(refractionV_Lens3L, normalV_Lens3R, Nn, Nn2)
            #ax.quiver(raySPoint_Lens3R[0],raySPoint_Lens3R[1],raySPoint_Lens3R[2],refractionV_Lens3R[0],refractionV_Lens3R[1],refractionV_Lens3R[2],linewidth=0.5)
            T = 0  # レンズ３とレンズ４の接着を考えた
            rayEPoint_Lens4L = raySPoint_Lens3R + T*refractionV_Lens3R
            VF.plotLineRed(raySPoint_Lens3R, rayEPoint_Lens4L)


            refractSPoint_Lens4L = rayEPoint_Lens4L  # 以下、レンズ４についての計算
            normalV_Lens4L = VF.decideNormalV_Lens4L(refractSPoint_Lens4L)
            #ax.quiver(refractSPoint_Lens4L[0],refractSPoint_Lens4L[1],refractSPoint_Lens4L[2],normalV_Lens4L[0],normalV_Lens4L[1],normalV_Lens4L[2],linewidth=0.5)
            refractionV_Lens4L = VF.decideRefractionVL(refractionV_Lens3R, normalV_Lens4L, Nn, Nn2)
            #ax.quiver(refractSPoint_Lens4L[0],refractSPoint_Lens4L[1],refractSPoint_Lens4L[2],refractionV_Lens4L[0],refractionV_Lens4L[1],refractionV_Lens4L[2],linewidth=0.5)
            T = VF.rayTraceDecideT_Lens4R(refractSPoint_Lens4L, refractionV_Lens4L)
            refractEPoint_Lens4L = refractSPoint_Lens4L + T*refractionV_Lens4L
            VF.plotLinePurple(refractSPoint_Lens4L,refractEPoint_Lens4L)  # 屈折光の描画
            raySPoint_Lens4R = refractEPoint_Lens4L
            normalV_Lens4R = VF.decideNormalV_Lens4R(raySPoint_Lens4R)
            #ax.quiver(raySPoint_Lens4R[0],raySPoint_Lens4R[1],raySPoint_Lens4R[2],normalV_Lens4R[0],normalV_Lens4R[1],normalV_Lens4R[2],linewidth=0.5)
            refractionV_Lens4R = VF.decideRefractionVR(refractionV_Lens4L, normalV_Lens4R, Nair, Nn2)
            #ax.quiver(raySPoint_Lens4R[0],raySPoint_Lens4R[1],raySPoint_Lens4R[2],refractionV_Lens4R[0],refractionV_Lens4R[1],refractionV_Lens4R[2],linewidth=0.5)
            T = 21  # 適当に置いた
            rayEPoint_Last = raySPoint_Lens4R + T*refractionV_Lens4R
            VF.plotLineBlue(raySPoint_Lens4R, rayEPoint_Last)

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
        ax.set_xlim(-LX+lens0V[0]+focusX, LX+lens0V[0]+focusX)
        ax.set_ylim(-LY+lens0V[1], LY+lens0V[1])
        ax.set_zlim(-LZ+lens0V[2], LZ+lens0V[2])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.title('Focus')


class ChromaticAberration:
    # 第０レンズと光線の描画
    def colorLens0(self):
        # 楕円面で球面レンズを再現する
        limitTheta = 2*np.pi  # theta生成数
        limitPhi = np.pi  # phi生成数
        theta = np.linspace(0, limitTheta, geneNum)
        phi = np.linspace(0, limitPhi, geneNum)
        Xs = Rx0 * np.outer(np.cos(theta), np.sin(phi)) + lens0V[0]
        Ys = Ry0 * np.outer(np.sin(theta), np.sin(phi)) + lens0V[1]
        Zs = Rz0 * np.outer(np.ones(np.size(theta)), np.cos(phi)) + lens0V[2]
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
            VF.plotLinePurple(raySPoint0, rayEPoint0)  # 入射光描画

            # 赤色光
            Nn = 1.5
            refractSPoint0 = rayEPoint0  # 入射光の終点を引き継ぐ
            normalV0 = VF.decideNormalV_Lens0(refractSPoint0)  # レンズの法線を求める
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],normalV0[0],normalV0[1],normalV0[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV0 = VF.decideRefractionVL(directionVector0, normalV0, Nair, Nn)
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],refractionV0[0],refractionV0[1],refractionV0[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens0R(refractSPoint0, refractionV0)
            refractEPoint0 = refractSPoint0 + T*refractionV0
            VF.plotLineRed(refractSPoint0,refractEPoint0)  # 屈折光の描画

            raySPoint1 = refractEPoint0  # 屈折光の終点を引き継ぐ
            normalV1 = VF.decideNormalV_Lens0(raySPoint1)  # レンズの法線を求める
            #ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],normalV1[0],normalV1[1],normalV1[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV1 = VF.decideRefractionVR(refractionV0, normalV1, Nair, Nn)
            #ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],refractionV1[0],refractionV1[1],refractionV1[2],linewidth=0.5)
            T = 12  # 衝突判定がないので適当に置いた
            rayEPoint1 = raySPoint1 + T*refractionV1  # 屈折光の終点
            VF.plotLineRed(raySPoint1,rayEPoint1)  # 屈折光の描画

            # 青色光、屈折率を赤色の1.04倍として計算
            Nn = Nn*1.04
            refractSPoint0 = rayEPoint0  # 入射光の終点を引き継ぐ
            normalV0 = VF.decideNormalV_Lens0(refractSPoint0)  # レンズの法線を求める
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],normalV0[0],normalV0[1],normalV0[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV0 = VF.decideRefractionVL(directionVector0, normalV0, Nair, Nn)
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],refractionV0[0],refractionV0[1],refractionV0[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens0R(refractSPoint0, refractionV0)
            refractEPoint0 = refractSPoint0 + T*refractionV0
            VF.plotLineBlue(refractSPoint0,refractEPoint0)  # 屈折光の描画

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

    # トリプレット（テッサー）と光線の描画。４枚目を低分散高屈折率の凸レンズにしている。
    def colorTriplet(self):
        def plotLensTriplet():
            # １枚目の凸レンズを再現する
            limitTheta = 2*np.pi  # theta生成数
            limitPhi = np.pi  # phi生成数
            theta = np.linspace(0, limitTheta, geneNum)
            phi = np.linspace(0, limitPhi, geneNum)
            Xs = Rx1 * np.outer(np.cos(theta), np.sin(phi))
            Xs = np.where(0<Xs, 0, Xs) + lens1V[0]
            Ys = Ry1 * np.outer(np.sin(theta), np.sin(phi)) + lens1V[1]
            Zs = Rz1 * np.outer(np.ones(np.size(theta)), np.cos(phi)) + lens1V[2]
            ax.plot_wireframe(Xs, Ys, Zs, linewidth=0.2)

            # ２枚目の凹レンズを再現する
            limitTheta = 2*np.pi  # theta生成数
            limitPhi = np.pi  # phi生成数
            theta = np.linspace(0, limitTheta, geneNum)
            phi = np.linspace(0, limitPhi, geneNum)
            Xs = Rx21 * np.outer(np.cos(theta), np.sin(phi))
            Ys = Ry21 * np.outer(np.sin(theta), np.sin(phi))
            Zs = Rz21 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Xs1 = Rx21 * np.outer(np.cos(theta), np.sin(phi))
            Xs2 = Rx22 * np.outer(np.cos(theta), np.sin(phi)) + 0.8
            Ys1 = Ry21 * np.outer(np.sin(theta), np.sin(phi))
            Ys2 = Ry22 * np.outer(np.sin(theta), np.sin(phi))
            Zs1 = Rz21 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Zs2 = Rz22 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Ys = np.where(Xs<0, Ys2, Ys1) + lens2V[1]
            Zs = np.where(Xs<0, Zs2, Zs1) + lens2V[2]
            Xs = np.where(Xs<0, Xs2, Xs1) + lens2V[0]
            ax.plot_wireframe(Xs, Ys, Zs, linewidth=0.2)

            # ３枚目の凹レンズを再現する
            limitTheta = 2*np.pi  # theta生成数
            limitPhi = np.pi  # phi生成数
            theta = np.linspace(0, limitTheta, geneNum)
            phi = np.linspace(0, limitPhi, geneNum)
            Xs = Rx31 * np.outer(np.cos(theta), np.sin(phi))
            Ys = Ry31 * np.outer(np.sin(theta), np.sin(phi))
            Zs = Rz31 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Xs1 = Rx31 * np.outer(np.cos(theta), np.sin(phi))
            Xs2 = Rx32 * np.outer(np.cos(theta), np.sin(phi)) + 1.1
            Ys1 = Ry31 * np.outer(np.sin(theta), np.sin(phi))
            Ys2 = Ry32 * np.outer(np.sin(theta), np.sin(phi))
            Zs1 = Rz31 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Zs2 = Rz32 * np.outer(np.ones(np.size(theta)), np.cos(phi))
            Ys = np.where(Xs<0, Ys2, Ys1) + lens3V[1]
            Zs = np.where(Xs<0, Zs2, Zs1) + lens3V[2]
            Xs = np.where(Xs<0, Xs2, Xs1) + lens3V[0] - 1.1
            ax.plot_wireframe(Xs, Ys, Zs, linewidth=0.2)

            # ４枚目のレンズを再現する
            limitTheta = 2*np.pi  # theta生成数
            limitPhi = np.pi  # phi生成数
            theta = np.linspace(0, limitTheta, geneNum)
            phi = np.linspace(0, limitPhi, geneNum)
            Xs = Rx4 * np.outer(np.cos(theta), np.sin(phi)) + lens4V[0]
            Ys = Ry4 * np.outer(np.sin(theta), np.sin(phi)) + lens4V[1]
            Zs = Rz4 * np.outer(np.ones(np.size(theta)), np.cos(phi)) + lens4V[2]
            ax.plot_wireframe(Xs, Ys, Zs, linewidth=0.2)

        plotLensTriplet()
        VF = VectorFunctions()  # インスタンス化

        # 始点を生成する
        pointsX = np.array([centerX]*25) + lens1V[0]
        pointsY, pointsZ = np.meshgrid(
            np.arange(-2+centerY, 3+centerY, 1),
            np.arange(-2+centerZ, 3+centerZ, 1))
        pointsY = pointsY.reshape(25)*rayDensity + lens1V[1]
        pointsZ = pointsZ.reshape(25)*rayDensity + lens1V[2]
        raySPoint0 = VF.makePoints(pointsX, pointsY, pointsZ, 25, 3)

        for i in raySPoint0:
            raySPoint0 = i
            directionVector0 = np.array([1, 0, 0])  # 入射光の方向ベクトルを設定
            T = VF.rayTraceDecideT_Lens1L(raySPoint0, directionVector0)  # 交点のための係数
            rayEPoint0 = raySPoint0 + T*directionVector0  # 入射光の終点
            VF.plotLinePurple(raySPoint0, rayEPoint0)  # 入射光描画

            # 赤色光
            Nn = 1.4
            Nn2 = 1.8
            refractSPoint0 = rayEPoint0  # 入射光の終点を引き継ぐ。以下レンズ１についての計算
            normalV_Lens1L = VF.decideNormalV_Lens1L(refractSPoint0)  # レンズの法線を求める
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],normalV_Lens1L[0],normalV_Lens1L[1],normalV_Lens1L[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens1L = VF.decideRefractionVL(directionVector0, normalV_Lens1L, Nair, Nn)
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],refractionV_Lens1L[0],refractionV_Lens1L[1],refractionV_Lens1L[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens1R(refractSPoint0, refractionV_Lens1L)
            refractEPoint0 = refractSPoint0 + T*refractionV_Lens1L
            VF.plotLineRed(refractSPoint0,refractEPoint0)  # 屈折光の描画
            raySPoint1 = refractEPoint0  # 屈折光の終点を引き継ぐ
            normalV1 = VF.decideNormalV_Lens1R(raySPoint1)  # レンズの法線を求める
            #ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],normalV1[0],normalV1[1],normalV1[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens1R = VF.decideRefractionVR(refractionV_Lens1L, normalV1, Nair, Nn)
            #ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],refractionV_Lens1R[0],refractionV_Lens1R[1],refractionV_Lens1R[2],linewidth=0.5)
            T = VF.rayTraceDecideT_Lens2L(raySPoint1, refractionV_Lens1R)
            rayEPoint1 = raySPoint1 + T*refractionV_Lens1R  # 空気中の屈折光の終点
            VF.plotLineRed(raySPoint1,rayEPoint1)  # 空気中の屈折光の描画

            refractSPoint_Lens2L = rayEPoint1  # 以下、レンズ２についての計算
            normalV_Lens2L = VF.decideNormalV_Lens2L(refractSPoint_Lens2L)  # レンズの法線を求める
            #ax.quiver(rayEPoint1[0],rayEPoint1[1],rayEPoint1[2],normalV_Lens2L[0],normalV_Lens2L[1],normalV_Lens2L[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens2L = VF.decideRefractionVL(refractionV_Lens1R, normalV_Lens2L, Nair, Nn)
            #ax.quiver(rayEPoint1[0],rayEPoint1[1],rayEPoint1[2],refractionV_Lens2L[0],refractionV_Lens2L[1],refractionV_Lens2L[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens2R(refractSPoint_Lens2L, refractionV_Lens2L)
            refractEPoint_Lens2L = refractSPoint_Lens2L + T*refractionV_Lens2L
            VF.plotLineRed(refractSPoint_Lens2L,refractEPoint_Lens2L)  # 屈折光の描画
            raySPoint_Lens2R = refractEPoint_Lens2L
            normalV_Lens2R = VF.decideNormalV_Lens2R(raySPoint_Lens2R)
            #ax.quiver(raySPoint_Lens2R[0],raySPoint_Lens2R[1],raySPoint_Lens2R[2],normalV_Lens2R[0],normalV_Lens2R[1],normalV_Lens2R[2],linewidth=0.5)
            refractionV_Lens2R = VF.decideRefractionVR(refractionV_Lens2L, normalV_Lens2R, Nair, Nn)
            #ax.quiver(raySPoint_Lens2R[0],raySPoint_Lens2R[1],raySPoint_Lens2R[2],refractionV_Lens2R[0],refractionV_Lens2R[1],refractionV_Lens2R[2],linewidth=0.5)
            T = VF.rayTraceDecideT_Lens3L(raySPoint_Lens2R, refractionV_Lens2R)
            rayEPoint_Lens3L = raySPoint_Lens2R + T*refractionV_Lens2R
            VF.plotLineRed(raySPoint_Lens2R, rayEPoint_Lens3L)

            refractSPoint_Lens3L = rayEPoint_Lens3L  # 以下、レンズ３についての計算
            normalV_Lens3L = VF.decideNormalV_Lens3L(refractSPoint_Lens3L)  # レンズの法線を求める
            #ax.quiver(refractSPoint_Lens3L[0],refractSPoint_Lens3L[1],refractSPoint_Lens3L[2],normalV_Lens3L[0],normalV_Lens3L[1],normalV_Lens3L[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens3L = VF.decideRefractionVL(refractionV_Lens2R, normalV_Lens3L, Nair, Nn)
            #ax.quiver(refractSPoint_Lens3L[0],refractSPoint_Lens3L[1],refractSPoint_Lens3L[2],refractionV_Lens3L[0],refractionV_Lens3L[1],refractionV_Lens3L[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens3R(refractSPoint_Lens3L, refractionV_Lens3L)
            refractEPoint_Lens3L = refractSPoint_Lens3L + T*refractionV_Lens3L
            VF.plotLineRed(refractSPoint_Lens3L,refractEPoint_Lens3L)  # 屈折光の描画
            raySPoint_Lens3R = refractEPoint_Lens3L
            normalV_Lens3R = VF.decideNormalV_Lens3R(raySPoint_Lens3R)
            #ax.quiver(raySPoint_Lens3R[0],raySPoint_Lens3R[1],raySPoint_Lens3R[2],normalV_Lens3R[0],normalV_Lens3R[1],normalV_Lens3R[2],linewidth=0.5)
            refractionV_Lens3R = VF.decideRefractionVR(refractionV_Lens3L, normalV_Lens3R, Nn, Nn2)
            #ax.quiver(raySPoint_Lens3R[0],raySPoint_Lens3R[1],raySPoint_Lens3R[2],refractionV_Lens3R[0],refractionV_Lens3R[1],refractionV_Lens3R[2],linewidth=0.5)
            T = 0  # レンズ３とレンズ４の接着を考えた
            rayEPoint_Lens4L = raySPoint_Lens3R + T*refractionV_Lens3R
            VF.plotLineRed(raySPoint_Lens3R, rayEPoint_Lens4L)


            refractSPoint_Lens4L = rayEPoint_Lens4L  # 以下、レンズ４についての計算
            normalV_Lens4L = VF.decideNormalV_Lens4L(refractSPoint_Lens4L)
            #ax.quiver(refractSPoint_Lens4L[0],refractSPoint_Lens4L[1],refractSPoint_Lens4L[2],normalV_Lens4L[0],normalV_Lens4L[1],normalV_Lens4L[2],linewidth=0.5)
            refractionV_Lens4L = VF.decideRefractionVL(refractionV_Lens3R, normalV_Lens4L, Nn, Nn2)
            #ax.quiver(refractSPoint_Lens4L[0],refractSPoint_Lens4L[1],refractSPoint_Lens4L[2],refractionV_Lens4L[0],refractionV_Lens4L[1],refractionV_Lens4L[2],linewidth=0.5)
            T = VF.rayTraceDecideT_Lens4R(refractSPoint_Lens4L, refractionV_Lens4L)
            refractEPoint_Lens4L = refractSPoint_Lens4L + T*refractionV_Lens4L
            VF.plotLineRed(refractSPoint_Lens4L,refractEPoint_Lens4L)  # 屈折光の描画
            raySPoint_Lens4R = refractEPoint_Lens4L
            normalV_Lens4R = VF.decideNormalV_Lens4R(raySPoint_Lens4R)
            #ax.quiver(raySPoint_Lens4R[0],raySPoint_Lens4R[1],raySPoint_Lens4R[2],normalV_Lens4R[0],normalV_Lens4R[1],normalV_Lens4R[2],linewidth=0.5)
            refractionV_Lens4R = VF.decideRefractionVR(refractionV_Lens4L, normalV_Lens4R, Nair, Nn2)
            #ax.quiver(raySPoint_Lens4R[0],raySPoint_Lens4R[1],raySPoint_Lens4R[2],refractionV_Lens4R[0],refractionV_Lens4R[1],refractionV_Lens4R[2],linewidth=0.5)
            T = 21  # 適当に置いた
            rayEPoint_Last = raySPoint_Lens4R + T*refractionV_Lens4R
            VF.plotLineRed(raySPoint_Lens4R, rayEPoint_Last)

            # 青色光
            Nn = Nn*1.04
            Nn2 = Nn2*1.01
            refractSPoint0 = rayEPoint0  # 入射光の終点を引き継ぐ。以下レンズ１についての計算
            normalV_Lens1L = VF.decideNormalV_Lens1L(refractSPoint0)  # レンズの法線を求める
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],normalV_Lens1L[0],normalV_Lens1L[1],normalV_Lens1L[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens1L = VF.decideRefractionVL(directionVector0, normalV_Lens1L, Nair, Nn)
            #ax.quiver(rayEPoint0[0],rayEPoint0[1],rayEPoint0[2],refractionV_Lens1L[0],refractionV_Lens1L[1],refractionV_Lens1L[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens1R(refractSPoint0, refractionV_Lens1L)
            refractEPoint0 = refractSPoint0 + T*refractionV_Lens1L
            VF.plotLineBlue(refractSPoint0,refractEPoint0)  # 屈折光の描画
            raySPoint1 = refractEPoint0  # 屈折光の終点を引き継ぐ
            normalV1 = VF.decideNormalV_Lens1R(raySPoint1)  # レンズの法線を求める
            #ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],normalV1[0],normalV1[1],normalV1[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens1R = VF.decideRefractionVR(refractionV_Lens1L, normalV1, Nair, Nn)
            #ax.quiver(raySPoint1[0],raySPoint1[1],raySPoint1[2],refractionV_Lens1R[0],refractionV_Lens1R[1],refractionV_Lens1R[2],linewidth=0.5)
            T = VF.rayTraceDecideT_Lens2L(raySPoint1, refractionV_Lens1R)
            rayEPoint1 = raySPoint1 + T*refractionV_Lens1R  # 空気中の屈折光の終点
            VF.plotLineBlue(raySPoint1,rayEPoint1)  # 空気中の屈折光の描画

            refractSPoint_Lens2L = rayEPoint1  # 以下、レンズ２についての計算
            normalV_Lens2L = VF.decideNormalV_Lens2L(refractSPoint_Lens2L)  # レンズの法線を求める
            #ax.quiver(rayEPoint1[0],rayEPoint1[1],rayEPoint1[2],normalV_Lens2L[0],normalV_Lens2L[1],normalV_Lens2L[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens2L = VF.decideRefractionVL(refractionV_Lens1R, normalV_Lens2L, Nair, Nn)
            #ax.quiver(rayEPoint1[0],rayEPoint1[1],rayEPoint1[2],refractionV_Lens2L[0],refractionV_Lens2L[1],refractionV_Lens2L[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens2R(refractSPoint_Lens2L, refractionV_Lens2L)
            refractEPoint_Lens2L = refractSPoint_Lens2L + T*refractionV_Lens2L
            VF.plotLineBlue(refractSPoint_Lens2L,refractEPoint_Lens2L)  # 屈折光の描画
            raySPoint_Lens2R = refractEPoint_Lens2L
            normalV_Lens2R = VF.decideNormalV_Lens2R(raySPoint_Lens2R)
            #ax.quiver(raySPoint_Lens2R[0],raySPoint_Lens2R[1],raySPoint_Lens2R[2],normalV_Lens2R[0],normalV_Lens2R[1],normalV_Lens2R[2],linewidth=0.5)
            refractionV_Lens2R = VF.decideRefractionVR(refractionV_Lens2L, normalV_Lens2R, Nair, Nn)
            #ax.quiver(raySPoint_Lens2R[0],raySPoint_Lens2R[1],raySPoint_Lens2R[2],refractionV_Lens2R[0],refractionV_Lens2R[1],refractionV_Lens2R[2],linewidth=0.5)
            T = VF.rayTraceDecideT_Lens3L(raySPoint_Lens2R, refractionV_Lens2R)
            rayEPoint_Lens3L = raySPoint_Lens2R + T*refractionV_Lens2R
            VF.plotLineBlue(raySPoint_Lens2R, rayEPoint_Lens3L)

            refractSPoint_Lens3L = rayEPoint_Lens3L  # 以下、レンズ３についての計算
            normalV_Lens3L = VF.decideNormalV_Lens3L(refractSPoint_Lens3L)  # レンズの法線を求める
            #ax.quiver(refractSPoint_Lens3L[0],refractSPoint_Lens3L[1],refractSPoint_Lens3L[2],normalV_Lens3L[0],normalV_Lens3L[1],normalV_Lens3L[2],linewidth=0.5)
            # 屈折光の方向ベクトルを求める
            refractionV_Lens3L = VF.decideRefractionVL(refractionV_Lens2R, normalV_Lens3L, Nair, Nn*0.88)
            #ax.quiver(refractSPoint_Lens3L[0],refractSPoint_Lens3L[1],refractSPoint_Lens3L[2],refractionV_Lens3L[0],refractionV_Lens3L[1],refractionV_Lens3L[2],linewidth=0.5)
            # 係数Tを求めて、屈折光の終点も求める
            T = VF.rayTraceDecideT_Lens3R(refractSPoint_Lens3L, refractionV_Lens3L)
            refractEPoint_Lens3L = refractSPoint_Lens3L + T*refractionV_Lens3L
            VF.plotLineBlue(refractSPoint_Lens3L,refractEPoint_Lens3L)  # 屈折光の描画
            raySPoint_Lens3R = refractEPoint_Lens3L
            normalV_Lens3R = VF.decideNormalV_Lens3R(raySPoint_Lens3R)
            #ax.quiver(raySPoint_Lens3R[0],raySPoint_Lens3R[1],raySPoint_Lens3R[2],normalV_Lens3R[0],normalV_Lens3R[1],normalV_Lens3R[2],linewidth=0.5)
            refractionV_Lens3R = VF.decideRefractionVR(refractionV_Lens3L, normalV_Lens3R, Nn*0.88, Nn2)
            #ax.quiver(raySPoint_Lens3R[0],raySPoint_Lens3R[1],raySPoint_Lens3R[2],refractionV_Lens3R[0],refractionV_Lens3R[1],refractionV_Lens3R[2],linewidth=0.5)
            T = 0  # レンズ３とレンズ４の接着を考えた
            rayEPoint_Lens4L = raySPoint_Lens3R + T*refractionV_Lens3R
            VF.plotLineBlue(raySPoint_Lens3R, rayEPoint_Lens4L)


            refractSPoint_Lens4L = rayEPoint_Lens4L  # 以下、レンズ４についての計算
            normalV_Lens4L = VF.decideNormalV_Lens4L(refractSPoint_Lens4L)
            #ax.quiver(refractSPoint_Lens4L[0],refractSPoint_Lens4L[1],refractSPoint_Lens4L[2],normalV_Lens4L[0],normalV_Lens4L[1],normalV_Lens4L[2],linewidth=0.5)
            refractionV_Lens4L = VF.decideRefractionVL(refractionV_Lens3R, normalV_Lens4L, Nn, Nn2)
            #ax.quiver(refractSPoint_Lens4L[0],refractSPoint_Lens4L[1],refractSPoint_Lens4L[2],refractionV_Lens4L[0],refractionV_Lens4L[1],refractionV_Lens4L[2],linewidth=0.5)
            T = VF.rayTraceDecideT_Lens4R(refractSPoint_Lens4L, refractionV_Lens4L)
            refractEPoint_Lens4L = refractSPoint_Lens4L + T*refractionV_Lens4L
            VF.plotLineBlue(refractSPoint_Lens4L,refractEPoint_Lens4L)  # 屈折光の描画
            raySPoint_Lens4R = refractEPoint_Lens4L
            normalV_Lens4R = VF.decideNormalV_Lens4R(raySPoint_Lens4R)
            #ax.quiver(raySPoint_Lens4R[0],raySPoint_Lens4R[1],raySPoint_Lens4R[2],normalV_Lens4R[0],normalV_Lens4R[1],normalV_Lens4R[2],linewidth=0.5)
            refractionV_Lens4R = VF.decideRefractionVR(refractionV_Lens4L, normalV_Lens4R, Nair, Nn2)
            #ax.quiver(raySPoint_Lens4R[0],raySPoint_Lens4R[1],raySPoint_Lens4R[2],refractionV_Lens4R[0],refractionV_Lens4R[1],refractionV_Lens4R[2],linewidth=0.5)
            T = 21  # 適当に置いた
            rayEPoint_Last = raySPoint_Lens4R + T*refractionV_Lens4R
            VF.plotLineBlue(raySPoint_Lens4R, rayEPoint_Last)

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
        ax.set_xlim(-LX/3+lens0V[0]+focusX, LX/3+lens0V[0]+focusX)
        ax.set_ylim(-LY/3+lens0V[1], LY/3+lens0V[1])
        ax.set_zlim(-LZ/3+lens0V[2], LZ/3+lens0V[2])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.title('Focus')



if __name__ == "__main__":
    start = time.time()

    CA = ChromaticAberration()  # インスタンス化
    fig = plt.figure(figsize=(16, 8))

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    CA.colorTriplet()  # レンズ付近を描画

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    CA.focusFocus(CA.colorTriplet)  # 焦点付近を描画

    print(time.time()-start)
    # グラフ描画
    plt.show()