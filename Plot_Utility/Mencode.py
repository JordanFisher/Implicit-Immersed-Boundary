import os

def Mencode(path, LastFrame, StepSize = 1, NameFormat = 'plot_%i.png', OutputFile = 'movie.avi'):
	"""Create a movie from image files.
	path is the directory path to the folder storing the images.
	LastFrame is the number of the last frame.
	StepSize is the difference in number between consecutive frames.
	NameFormat gives the filename of each image, %i is a wildcard for the frame number.
	OutputFile is the filename of the .avi created"""
	
	name = NameFormat + '\n'
	
	# Set cwd (current working directory)
	os.chdir(path)
	
	# Create list of images to use, saving the list to "list.txt"
	outfile = open(os.path.join(os.getcwd(), "list.txt"), 'w')

	for i in range(StepSize, LastFrame, StepSize):
		outfile.write(name % i)
		
	outfile.close()

	
	# Setup mencoder command
	mencoder = 'C:\\MPlayer\\mencoder' # The location of mencoder
	files = 'mf://@list.txt'           # Tell mencoder to use "list.txt"

	command = '%s %s -mf fps=120 -o %s -ovc xvid -xvidencopts bitrate=3000' % (mencoder, files, OutputFile)

	print command 
	print

	# Use the command to create the movie
	os.system(command)