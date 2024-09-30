import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

def InitGaussFilter(this):
	this['TapCount'] = this['symbol span'] * this['samples per symbol']
	# generate normalized time
	this['Time'] = np.arange(0, this['symbol span'], 1 / this['samples per symbol']) - ((this['symbol span'] - (1 / this['samples per symbol'])) / 2)
	this['SymbolTicks'] = np.arange(-this['symbol span'] / 2, this['symbol span'] / 2, 1)
	this['Taps'] = np.zeros(this['TapCount'])

	alpha = pow(np.log(2) / 2, 0.5) / this['BT']
	print('Gaussian alpha: ', round(alpha,2))
	index = 0
	for time in this['Time']:
		numerator = np.exp(-pow(np.pi,2) * pow(time, 2) / pow(alpha, 2))
		try:
			this['Taps'][index] = numerator
		except:
			pass
		index += 1

	
	this['Taps'] = np.convolve(this['Taps'], np.ones(this['samples per symbol']), 'same')
	this['Taps'] = this['Taps'] / max(this['Taps'])
	this['segment count'] = np.power(2,filter['symbol span'])
	return this
	
def GenerateSegment(this, data):
	data = int(data)
	data_stream = np.zeros((this['symbol span']) * this['samples per symbol'])
	stream_index = 0
	# decompose the input data into binary symbols in a list
	# Create a tap for the MSB of the symbol span
	data_tap = np.power(2, this['symbol span'] - 1)
	symbol_stream = []
	for symbol_index in range(this['symbol span']):
		if data & data_tap:
			symbol_stream.append(-1)
		else:
			symbol_stream.append(1)
		data <<= 1
	#print(symbol_stream)
	
	for symbol_index in range(this['symbol span']):
		#print(symbol_index)
		data_stream[stream_index] = symbol_stream[symbol_index]
		stream_index += this['samples per symbol']
	response = np.convolve(data_stream, this['Taps'], 'full')
	first_valid_index = np.floor((len(this['Taps']) - 1) / 2)
	first_valid_index += np.floor((this['symbol span'] * this['samples per symbol'] - 1) / 2) - this['samples per symbol']
	first_valid_index = int(first_valid_index)
	last_valid_index = first_valid_index + this['samples per symbol']
	last_valid_index = int(last_valid_index)
	#print(f"first_valid_index {first_valid_index}, last valid index {last_valid_index}")
	response = response[first_valid_index:last_valid_index]
	return(response)
	
def GenInt16ArrayC(name, array, column_width):
	result = '\n'
	result += f'const __prog__ int16_t __attribute__((space(prog))) {name}[{len(array)}] = '
	result += '{ '
	y = len(array)
	for x in range(y):
		if x % column_width == 0:
			result += ' \\\n     '
		result += f' {int(np.rint(array[x]))}'
		if x < (y-1):
			result += ','
	result += ' };'
	return result

if len(sys.argv) != 6:
		print("Not enough arguments. Usage: python3 gfsk-gen.py <bandwidth-time> <span> <samples per symbol> <center freq> <deviation>")
		sys.exit(-1)

filter = {}
filter['BT'] = float(sys.argv[1])
filter['symbol span'] = int(sys.argv[2])
filter['samples per symbol'] = int(sys.argv[3])
filter['center freq'] = float(sys.argv[4])
filter['deviation'] = float(sys.argv[5])

filter = InitGaussFilter(filter)

plt.figure()
plt.plot(filter['Time'],filter['Taps'])
plt.xticks(filter['SymbolTicks'])
plt.title(f'Gaussian Pulse Filter, BT: {filter["BT"]}')
plt.grid()
plt.show()

# now generate segments of the baseband waveform
waveform_table = []
print(f"segment count: {filter['segment count']}")
for index in range(filter['segment count']):
	segment = GenerateSegment(filter, index)
	waveform_table.extend(segment)

# 
print(f"length of segment: {len(waveform_table)}")

for index in range(len(waveform_table)):
	waveform_table[index] *= filter['deviation']
	waveform_table[index] += filter['center freq']
	waveform_table[index] = int(round(waveform_table[index])) 

plt.figure()
half_point = len(waveform_table) / 2
half_point = int(np.floor(half_point) - 1)
ymin = filter['center freq'] - (2 * filter['deviation'])
ymax = filter['center freq'] + (2 * filter['deviation'])
plt.plot(waveform_table[0:half_point])
plt.plot(waveform_table[half_point+1:])
plt.ylim(ymin, ymax)
plt.show()


#generate a new director for the reports
run_number = 0
print('trying to make a new directory')
while True:
	run_number = run_number + 1
	dirname = f'./run{run_number}/'
	try:
		os.mkdir(dirname)
	except:
		print(dirname + ' exists')
		continue
	print(f'made new directory {dirname}')
	break

# Generate and save report file
report_file_name = f'run{run_number}_report.txt'
try:
	report_file = open(dirname + report_file_name, 'w+')
except:
	print('Unable to create report file.')
with report_file:
	report_file.write('# Command line: ')
	for argument in sys.argv:
		report_file.write(f'{argument} ')

	report_file.write('\n\n# Gaussian Pulse Filter\n')
	report_file.write('\n')
	report_file.write(GenInt16ArrayC(f'GaussFilter', filter['Taps'] * 32768, filter['samples per symbol']))
	report_file.write('\n')
	report_file.write(GenInt16ArrayC(f'Waveform table', waveform_table, filter['samples per symbol']))
	report_file.write('\n')
	report_file.close()
	print(f'wrote {dirname + report_file_name}')
