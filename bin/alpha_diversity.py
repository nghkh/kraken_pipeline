#!/usr/bin/env python
import os, sys, argparse
import math

def shannons_alpha(p):
	# shannons = -sum(pi ln(pi))
	h = []
	for i in p:
		if i > 0:  # Avoid math domain error with log(0)
			h.append(i * math.log(i))
	print("Shannon's diversity: %s" %(-1 *sum(h)))
	return (-1 *sum(h))

def berger_parkers_alpha(p):
	# bp is nmax/N which is equal to the max pi == max(p)
	if not p:  # Check if p is empty
		print("Berger-parker's diversity: NA (no valid data)")
		return 0
	print("Berger-parker's diversity: %s" %max(p))
	return max(p)

def simpsons_alpha(D):
	# simpsons index of diversity = 1 - D
	# D = (sum ni(ni-1))/ N*(N-1)
	print("Simpson's index of diversity: %s" %(1-D))
	return 1-D

def inverse_simpsons_alpha(D):
	# simsons inverse = 1/D
	if D == 0:  # Avoid division by zero
		print("Simpson's Reciprocal Index: NA (D=0)")
		return float('inf')
	print("Simpson's Reciprocal Index: %s" %(1/D))
	return 1/D

def fishers_alpha():	
	try:
		global np
		import numpy as np
		from scipy.optimize import fsolve
		
		# Check if we have enough data for meaningful calculation
		if S_f <= 1 or N_f <= S_f:
			print("Fisher's index: NA (insufficient data)")
			return [0]
			
		fish = fsolve(eqn_output, 1)
		
		# NOTE:
		# if ratio of N/S > 20 then x > 0.99 (Poole,1974)
		# x is almost always > 0.9 and never > 1.0 i.e. ~ 0.9 < x < 1.0
		print("Fisher's index: %s" %fish[0])
		return fish
	except ImportError:
		print("Fisher's index: NA (numpy or scipy not available)")
		return [0]
	except:
		print("Fisher's index: NA (calculation error)")
		return [0]

def eqn_output(a):
	return a * np.log(1+N_f/a) - S_f

def parse_bracken_file(filename):
	"""Parse a Bracken file format"""
	try:
		f = open(filename)
		header = f.readline()  # Skip header
		n = []
		
		# Bracken output format has abundance info in the "new_est_reads" column
		# Which is typically column 5 (0-indexed)
		for line in f:
			cols = line.strip().split('\t')
			if len(cols) < 6:  # Standard Bracken output should have at least 6 columns
				continue
				
			try:
				# Get abundance from the "new_est_reads" column 
				ind_abund = float(cols[5])
				if ind_abund > 0:  # Only include non-zero values
					n.append(ind_abund)
			except (ValueError, IndexError):
				continue
				
		f.close()
		
		# Check if we found any data
		if not n:
			print(f"Warning: No abundance values found in {filename}")
			with open(filename, 'r') as debug_f:
				print(f"File content preview: {debug_f.read(500)}")
				
		return n
	except Exception as e:
		print(f"Error parsing Bracken file: {str(e)}")
		import traceback
		traceback.print_exc()
		return []

# Main method
def main():
	# get arguments
	parser = argparse.ArgumentParser(description='Calculate alpha diversity from Bracken abundance estimation file.')
	parser.add_argument('-f','--filename', dest='filename', required=True, 
		help='Bracken output file')
	parser.add_argument('-a','--alpha', dest='value', required=True, type=str,  
		help='Single letter alpha diversity type (Sh, BP, Si, ISi, F)')
	parser.add_argument('--type', dest='filetype', default='bracken', 
		choices=['bracken'],
		help='Type of input file: bracken. Default = bracken')
	parser.add_argument('--level', dest='level', default='S',
		help='Taxonomic level from which the Bracken file was created')
	args = parser.parse_args()

	# Check if file exists
	if not os.path.isfile(args.filename):
		print(f"Error: File {args.filename} not found")
		sys.exit(1)
	
	# Parse Bracken file
	n = parse_bracken_file(args.filename)
	
	# Check if we have data
	if not n:
		print(f"Error: No valid abundance data found in {args.filename}")
		sys.exit(1)
	
	# calculations
	N = sum(n)  # total number of individuals
	if N == 0:
		print(f"Error: Total abundance is zero in {args.filename}")
		sys.exit(1)
		
	S = len(n)  # total number of species
	
	# calculate all the pi's
	p = []  # store pi's 
	D = 0
	for i in n:  # go through each species
		if i != 0:  # there should not be any zeros
			pi = i/N  # pi is the ni/N
			p.append(pi)
			D += i*(i-1)
	
	if N > 1:
		D = D/(N*(N-1))
	else:
		D = 0
	
	# find the indicated alpha
	if args.value == 'Sh':  # calculate shannon's diversity
		result = shannons_alpha(p)
	elif args.value == 'BP':  # calculate berger-parker's dominance index
		result = berger_parkers_alpha(p)
	elif args.value == 'Si':  # calculate Simpson's alpha 
		result = simpsons_alpha(D)
	elif args.value == 'ISi':  # calculate Inverse Simpson's alpha 
		result = inverse_simpsons_alpha(D)
	elif args.value == 'F':  # calculate fisher's alpha
		print("Fisher's alpha...loading")
		global N_f
		N_f = N
		global S_f
		S_f = S
		result = fishers_alpha()
	else:
		print("Not a supported alpha")
		sys.exit(1)

if __name__ == "__main__":
    main()

