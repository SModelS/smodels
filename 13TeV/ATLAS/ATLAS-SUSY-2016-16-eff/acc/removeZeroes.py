from os import listdir, getcwd
from os.path import isfile, join

path = getcwd()
effmaps = [f for f in listdir(path) if isfile(join(path, f)) and f[:6] == 'EffMap']

for em in effmaps:

	print('file: {}'.format(em))

	with open(em, 'r') as f:
		lines = f.readlines()

	newlines = []
	first 	 = False
	last  	 = ''

	for line in lines:

		if line[0] == '#':
			newlines.append(line)
		else:
			split = line.split(' ')

			if float(split[2][:-1]) != 0.0:
				newlines.append(line)
			elif first == False:
				newlines.append(line)
				first = True
			else:
				last = line

	newlines.append(last)
	#for line in newlines: print(line)

	with open(em, 'w') as f:
		for line in newlines:
			f.write(line)

	print('  ..done!')
