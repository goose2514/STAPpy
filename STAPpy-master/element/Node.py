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
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

class CNode(object):
	NDF = 3

	def __init__(self, x=0.0, y=0.0, z=0.0):
		super().__init__()
		self.XYZ = np.zeros(CNode.NDF)
		self.XYZ[0] = x
		self.XYZ[1] = y
		self.XYZ[2] = z

		self.bcode = np.zeros(CNode.NDF, dtype=int)
		self.NodeNumber = 0

	def Read(self, input_file, check_np):
		line = input_file.readline().split()

		N = int(line[0])
		if N != check_np + 1:
			error_info = "\n*** Error *** Nodes must be inputted in order !" \
						 "\n   Expected node number : {}" \
						 "\n   Provided node number : {}".format(check_np+1, N)
			raise ValueError(error_info)

		self.NodeNumber = N

		self.bcode[0] = int(line[1])
		self.bcode[1] = int(line[2])
		self.bcode[2] = int(line[3])
		self.XYZ[0] = np.double(line[4])
		self.XYZ[1] = np.double(line[5])
		self.XYZ[2] = np.double(line[6])

	def Write(self, output_file):
		node_info = "%9d%5d%5d%5d%18.6e%15.6e%15.6e\n"%(
			self.NodeNumber, self.bcode[0], self.bcode[1], self.bcode[2],
			self.XYZ[0], self.XYZ[1], self.XYZ[2])
		print(node_info, end='')
		output_file.write(node_info)

	def WriteEquationNo(self, output_file):
		equation_info = "%9d       "%self.NodeNumber

		for dof in range(CNode.NDF):
			equation_info += "%5d"%self.bcode[dof]

		equation_info += '\n'
		print(equation_info, end='')
		output_file.write(equation_info)

	def WriteNodalDisplacement(self, output_file, displacement):
		displacement_info = "%5d        "%self.NodeNumber

		for dof in range(CNode.NDF):
			if self.bcode[dof] == 0:
				displacement_info += "%18.6e"%0.0
			else:
				displacement_info += "%18.6e"%displacement[self.bcode[dof] - 1]

		displacement_info += '\n'
		print(displacement_info, end='')
		output_file.write(displacement_info)