from numpy import array, zeros, float64

import GetDecomposition
import Panel as _Panel
		
class TreeCode:
	def __init__(self, fluid, Terms = 10, MinPointsToHaveChildren = 10, MinPanelWidth = 2, Bounds = [[0,0,0],[1,1,1]]):
		self.Terms, self.MinPointsToHaveChildren, self.MinPanelWidth, self.Bounds = Terms, MinPointsToHaveChildren, MinPanelWidth, Bounds
	
		self.Bounds = array(self.Bounds, float64)
	
		# Initialize the lookup tables needed by the tree code
		
		# Get the fluid Greens function lookup table
		self.UField = fluid.GetLookup()

		# Initialize the far field decomposition lookup table
		self.DecompGroup = GetDecomposition.GetDecomposition(fluid, Terms = Terms, MinPanelWidth = MinPanelWidth)
		
	def Initialize(self, X):
		Terms, MinPointsToHaveChildren, MinPanelWidth, Bounds = self.Terms, self.MinPointsToHaveChildren, self.MinPanelWidth, self.Bounds
		DecompGroup, UField = self.DecompGroup, self.UField
		
		# Initialize the tree structure
		self.TreeRoot = _Panel.Panel(MinPointsToHaveChildren, DecompGroup.Depth, X, Bounds[0], Bounds[1], 0)
		
		self.TreeRoot.CalcUV(Terms, DecompGroup)
		self.TreeRoot.PreEval(UField, DecompGroup, X)
	

	def FastFluidSolve(self, fiber):
		"""Perform a fast fluid solve utilizing the treecode.
		
		The treecode evaluation includes the spreading of the fiber force and the interpolation back to the fiber."""
		
		self.TreeRoot.CalcSeriesTerms(fiber.F)
		self.TreeRoot.FastEval(fiber.LaggedX, fiber.F, fiber.U)
		
		fiber.U *= fiber.hb
		
		return fiber.U