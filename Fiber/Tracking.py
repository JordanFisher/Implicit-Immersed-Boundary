import pylab as m
from numpy import dot, array
from IBUtility import EnsureNumpyArray, cross
from Fiber import Fiber

class Tracking(Fiber):
	"""A collection of non-interacting fiber points used to track fluid flow."""
	
	def __init__(self, fluid, N=0, Dim=2, Mask=None, Groups=3):
		"""N is the number of fiber points to create along one length of the fluid domain."""

		self.fluid = fluid
		
		Fiber.__init__(self, N**2, Dim = Dim)

		h = fluid.dims / N

		self.GroupStartIndices = GroupStartIndices = [0]
		
		GroupWidth = N / Groups
		if GroupWidth == 0:
			GroupWidth = Groups = 1
		NumOfBigGroups = N % Groups # The remainder is the number of groups that must be one column wider
			
		index = 0
		EndOfGroup = GroupWidth
		CurGroup = 0
		for j in range(N):
			# Check to see if we need to start a new group
			if j == EndOfGroup:
				CurGroup += 1
				GroupStartIndices.append(index)
				
				# Calculate the ending index of the new group
				if CurGroup < Groups - NumOfBigGroups:
					EndOfGroup += GroupWidth
				else:
					EndOfGroup += GroupWidth + 1				
			
			# Loop through all the points in this column
			for k in range(N):
				# Calculate the physical position of the current point
				x, y = (j + .5) * h[0], (k + .5) * h[1]
				
				# Check to see if the point is masked
				if Mask != None:
					n, m = int(x / fluid.h[0]), int(y / fluid.h[1])
					if Mask[n,m] != 0:
						continue				
				
				self.X[index,:] = [x, y]
				
				index += 1
				
		GroupStartIndices.append(index)
		self.Groups = Groups
				
		hold = self.X
		Fiber.__init__(self, index, Dim = Dim)
		self.X = hold[:self.Nb]
		
	def Update(self):
		X = self.X
		
		self.FromGrid(self.fluid)
		X += self.fluid.dt * self.U
		
		for i in range(self.Dim):
			X[...,i] = X[...,i] % self.fluid.dims[i]
	
	def Plot(self):
		X = self.X
		
		Indices = self.GroupStartIndices
		
		c = ['g', 'b', 'r']
		
		for i in range(self.Groups):	
			i1, i2 = Indices[i], Indices[i+1]
			m.scatter(X[i1:i2,0], X[i1:i2,1], marker = 's', c = c[i % len(c)])

		self.fluid.SetPlotExtent()
				