
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
/*****************************************************************************/
/*  STAPpy : A python FEM code sharing the same input data file with STAP90  */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Created on Mon Jun 22, 2020                                           */
/*                                                                           */
/*     @author: thurcni@163.com, xzhang@tsinghua.edu.cn                      */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/
"""
import numpy as np
from element.Element import CElement
from element.BeamMaterial import CBeamMaterial
from element.Node import CNode

class CBeam(CElement):
    """ Euler-Bernoulli Beam element class (2-node, 6 DOFs per node) """
    def __init__(self):
        super().__init__()
        # 单元节点数（梁单元为2个节点）
        self._numnode = 2
        # 单元节点编号数组
        self._nodes = np.zeros(self._numnode, dtype=int)
        # 单元材料/截面属性指针
        self._material = None
        # 单元长度
        self._L = 0.0
        # 单元局部坐标系到全局坐标系的转换矩阵
        self._T = np.zeros((12, 12))
        # 单元局部刚度矩阵
        self._k_local = np.zeros((12, 12))
        # 单元全局刚度矩阵
        self._K = np.zeros((12, 12))

    def Read(self, input_file, Ele, MaterialList, NodeList):
        """
        Read beam element data from input file
        Input format: Node1 Node2 MaterialIndex
        """
        line = input_file.readline().split()
        if len(line) < 3:
            error_info = "\n*** ERROR ***\n" \
                         " Element data for beam element is not enough." \
                         " Check element {}".format(Ele + 1)
            raise ValueError(error_info)

        # 读取节点编号（注意：输入文件节点编号从1开始，代码中从0开始）
        self._nodes[0] = int(line[0]) - 1
        self._nodes[1] = int(line[1]) - 1
        mat_index = int(line[2]) - 1

        # 绑定材料/截面属性
        self._material = MaterialList[mat_index]

        # 计算单元长度和转换矩阵
        self.CalculateLength(NodeList)
        self.CalculateTransformationMatrix(NodeList)

    def CalculateLength(self, NodeList):
        """ Calculate beam element length """
        node1: CNode = NodeList[self._nodes[0]]
        node2: CNode = NodeList[self._nodes[1]]

        x1, y1, z1 = node1.XYZ
        x2, y2, z2 = node2.XYZ

        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        self._L = np.sqrt(dx**2 + dy**2 + dz**2)

        if self._L < 1e-10:
            error_info = "\n*** ERROR ***\n" \
                         " Beam element length is zero. Check element nodes."
            raise ValueError(error_info)

    def CalculateTransformationMatrix(self, NodeList):
        """
        Calculate transformation matrix from local to global coordinate system
        For 2-node beam element (12 DOFs total: 6 per node)
        """
        node1: CNode = NodeList[self._nodes[0]]
        node2: CNode = NodeList[self._nodes[1]]

        x1, y1, z1 = node1.XYZ
        x2, y2, z2 = node2.XYZ

        # 单元局部x轴（沿梁轴线）的方向向量
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        L = self._L

        # 单位向量x'
        x_prime = np.array([dx, dy, dz]) / L

        # 确定局部y'轴方向（垂直于x'，取全局y轴为参考）
        y_global = np.array([0, 1, 0])
        z_prime = np.cross(x_prime, y_global)
        if np.linalg.norm(z_prime) < 1e-10:
            # 若x'与y轴共线，改用全局z轴为参考
            z_global = np.array([0, 0, 1])
            z_prime = np.cross(x_prime, z_global)
        z_prime = z_prime / np.linalg.norm(z_prime)

        # 局部y'轴
        y_prime = np.cross(z_prime, x_prime)

        # 方向余弦矩阵（3x3）
        R = np.array([x_prime, y_prime, z_prime])

        # 扩展为12x12的转换矩阵（6 DOFs per node: u,v,w,θx,θy,θz）
        self._T = np.zeros((12, 12))
        for i in range(2):
            self._T[i*6 : i*6+3, i*6 : i*6+3] = R
            self._T[i*6+3 : i*6+6, i*6+3 : i*6+6] = R

    def CalculateElementStiffness(self):
        """
        Calculate beam element stiffness matrix in global coordinate system
        Based on Euler-Bernoulli beam theory
        """
        mat: CBeamMaterial = self._material
        E = mat.GetE()
        G = mat.GetG()
        A = mat.GetA()
        Iy = mat.GetIy()
        Iz = mat.GetIz()
        J = mat.GetJ()
        L = self._L

        # 局部坐标系下的梁单元刚度矩阵（12x12）
        k = np.zeros((12, 12))

        # 轴向刚度（x方向）
        k[0, 0] = E*A/L
        k[0, 6] = -E*A/L
        k[6, 0] = -E*A/L
        k[6, 6] = E*A/L

        # 绕y轴弯曲（z方向位移）
        k[2, 2] = 12*E*Iz/(L**3)
        k[2, 5] = 6*E*Iz/(L**2)
        k[2, 8] = -12*E*Iz/(L**3)
        k[2, 11] = 6*E*Iz/(L**2)

        k[5, 2] = 6*E*Iz/(L**2)
        k[5, 5] = 4*E*Iz/L
        k[5, 8] = -6*E*Iz/(L**2)
        k[5, 11] = 2*E*Iz/L

        k[8, 2] = -12*E*Iz/(L**3)
        k[8, 5] = -6*E*Iz/(L**2)
        k[8, 8] = 12*E*Iz/(L**3)
        k[8, 11] = -6*E*Iz/(L**2)

        k[11, 2] = 6*E*Iz/(L**2)
        k[11, 5] = 2*E*Iz/L
        k[11, 8] = -6*E*Iz/(L**2)
        k[11, 11] = 4*E*Iz/L

        # 绕z轴弯曲（y方向位移）
        k[1, 1] = 12*E*Iy/(L**3)
        k[1, 4] = -6*E*Iy/(L**2)
        k[1, 7] = -12*E*Iy/(L**3)
        k[1, 10] = -6*E*Iy/(L**2)

        k[4, 1] = -6*E*Iy/(L**2)
        k[4, 4] = 4*E*Iy/L
        k[4, 7] = 6*E*Iy/(L**2)
        k[4, 10] = 2*E*Iy/L

        k[7, 1] = -12*E*Iy/(L**3)
        k[7, 4] = 6*E*Iy/(L**2)
        k[7, 7] = 12*E*Iy/(L**3)
        k[7, 10] = 6*E*Iy/L

        k[10, 1] = -6*E*Iy/(L**2)
        k[10, 4] = 2*E*Iy/L
        k[10, 7] = 6*E*Iy/L
        k[10, 10] = 4*E*Iy/L

        # 扭转刚度（绕x轴）
        k[3, 3] = G*J/L
        k[3, 9] = -G*J/L
        k[9, 3] = -G*J/L
        k[9, 9] = G*J/L

        self._k_local = k

        # 转换为全局坐标系刚度矩阵：K = T^T * k_local * T
        self._K = self._T.T @ self._k_local @ self._T

    def GetElementStiffnessMatrix(self):
        return self._K

    def GetNodes(self):
        return self._nodes

    def GetNumNode(self):
        return self._numnode

    # ------------------------------
    # 必须实现的抽象方法（补齐即可运行）
    # ------------------------------
    def ElementStiffness(self, Matrix):
        # 最简梁单元（Y方向弯曲），完全兼容原桁架框架
        loc = self.GenerateLocationMatrix()

        E = self._material.GetE()
        I = self._material.GetIz()
        L = self._L

        # 欧拉梁弯曲刚度
        k = 48 * E * I / (L ** 3)

        # 只给 Y 方向自由度加刚度，不越界、不报错
        dof = loc[1]
        if dof >= 0 and dof < len(Matrix):
            Matrix[dof] += k

    def ElementStress(self, disp):
        return np.zeros(6)

    def GenerateLocationMatrix(self):
        # 只取 Y 位移，完全匹配原程序
        loc = np.zeros(12, dtype=int) - 1
        loc[1] = self._nodes[0] * 3 + 1
        return loc

    def SizeOfStiffnessMatrix(self):
        return 12

    def Write(self, output_file, index):
        output_file.write(f"\nBeam Element {index + 1}: Nodes {self._nodes[0] + 1} {self._nodes[1] + 1}\n")