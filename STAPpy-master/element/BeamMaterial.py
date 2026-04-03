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
from element.Material import CMaterial

class CBeamMaterial(CMaterial):
    """ Beam element material and section property class """
    def __init__(self):
        super().__init__()
        # 弹性模量
        self._E = 0.0
        # 剪切模量
        self._G = 0.0
        # 截面面积
        self._A = 0.0
        # 绕y轴惯性矩
        self._Iy = 0.0
        # 绕z轴惯性矩
        self._Iz = 0.0
        # 扭转惯性矩
        self._J = 0.0

    def Read(self, input_file, mset):
        """
        Read beam material and section property from input file
        Input format: E G A Iy Iz J
        """
        line = input_file.readline().split()
        if len(line) < 6:
            error_info = "\n*** ERROR ***\n" \
                         " Material/section property for beam element is not enough." \
                         " Check line {}".format(mset + 1)
            raise ValueError(error_info)

        self._E = float(line[0])
        self._G = float(line[1])
        self._A = float(line[2])
        self._Iy = float(line[3])
        self._Iz = float(line[4])
        self._J = float(line[5])

    def GetE(self):
        return self._E

    def GetG(self):
        return self._G

    def GetA(self):
        return self._A

    def GetIy(self):
        return self._Iy

    def GetIz(self):
        return self._Iz

    def GetJ(self):
        return self._J

    def Write(self, output_file, mset):
        """
        Write beam material data to output file
        """
        output_file.write(f"\nMaterial {mset + 1}:\n")
        output_file.write(f"  E  = {self._E:.6e}\n")
        output_file.write(f"  G  = {self._G:.6e}\n")
        output_file.write(f"  A  = {self._A:.6e}\n")
        output_file.write(f"  Iy = {self._Iy:.6e}\n")
        output_file.write(f"  Iz = {self._Iz:.6e}\n")
        output_file.write(f"  J  = {self._J:.6e}\n")