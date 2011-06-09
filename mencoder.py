import os

#path = 'C:\\Dissertation\\Simulation_Data\\Peristaltic\\Explicit2D'
path = 'C:\\Dissertation\\Simulation_Data\\Peristaltic\\Explicit2D_NoKill_RK2_W5_N128'

name = "plot_%i.png\n"
output = 'movie.avi'

end, step = 262500, 500


# Set cwd
os.chdir(path)

# Create list of images to use
outfile = open(os.path.join(os.getcwd(), "list.txt"), 'w')

for i in range(step, end, step):
    outfile.write(name % i)
	
outfile.close()


# Create movie using mencoder
mencoder = 'C:\\MPlayer\\mencoder'
files = 'mf://@list.txt'

command = '%s %s -mf fps=120 -o %s -ovc xvid -xvidencopts bitrate=3000' % (mencoder, files, output)

print command 
print

##(dummy, stdout_and_stderr) = os.popen4(command)
##print stdout_and_stderr.read()

os.system(command)

##C:\MPlayer\mencoder "mf://*.png" -mf fps=30 -o output.avi -ovc xvid -xvidencopts bitrate=3000
