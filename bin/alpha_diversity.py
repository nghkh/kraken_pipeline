#!/usr/bin/env python
import os, sys, argparse
import math

def shannons_alpha(p):
	# shannons = -sum(pi ln(pi))
	h = []
	for i in p:
		h.append(i * math.log(i))
	print("Shannon's diversity: %s" %(-1 *sum(h)))
	return (-1 *sum(h))

def berger_parkers_alpha(p):
	# bp is nmax/N which is equal to the max pi == max(p)
	print("Berger-parker's diversity: %s" %max(p))
	return max(p)

def simpsons_alpha(D):
	# simpsons index of diversity = 1 - D
	# D = (sum ni(ni-1))/ N*(N-1)
	print("Simpson's index of diversity: %s" %(1-D))
	return 1-D

def inverse_simpsons_alpha(D):
	# simsons inverse = 1/D
	print("Simpson's Reciprocal Index: %s" %(1/D))
	return 1/D

def fishers_alpha():	
	global np
	import numpy as np
	from scipy.optimize import fsolve
	
	fish = fsolve(eqn_output,1)
	
	# NOTE:
	# if ratio of N/S > 20 then x > 0.99 (Poole,1974)
	# x is almost always > 0.9 and never > 1.0 i.e. ~ 0.9 < x < 1.0
	print("Fisher's index: %s" %fish[0])
	return fish

def eqn_output(a):
	return a * np.log(1+N_f/a) - S_f

# Main method
def main():
	# get arguments
	parser = argparse.ArgumentParser(description='Calculate alpha diversity metrics.')
	parser.add_argument('-f','--filename', dest='filename', help='Input file with abundance estimates (Bracken, Kraken, etc.)')
	parser.add_argument('-a','--alpha', dest='value', default='Sh', type=str, 
                       help='Type of alpha diversity to calculate: Sh (Shannon), BP (Berger-Parker), Si (Simpson), ISi (Inverse Simpson), F (Fisher). Default: Sh')
	parser.add_argument('--type', dest='filetype', default='bracken', choices=['bracken', 'kreport', 'kreport2', 'krona'],
                       help='Type of input file: bracken, kreport, kreport2, krona. Default: bracken')
	parser.add_argument('--level', dest='level', default='all', choices=['all', 'D', 'P', 'C', 'O', 'F', 'G', 'S'],
                       help='Taxonomic level to analyze. Default: all')
	args = parser.parse_args()

	f = open(args.filename)
	# Skip header line if present
	header = f.readline()
	
	n = []  # Store abundances
	tax_col = 0   # Column containing taxonomy info
	count_col = 1  # Column containing count info
	level_col = -1  # Column containing taxonomy level info
	
	# Set columns based on file type
	if args.filetype == 'bracken':
		count_col = 5    # Abundance column in Bracken
		tax_col = 0      # Taxonomy column in Bracken
		level_col = 2    # Level column in Bracken
	elif args.filetype in ['kreport', 'kreport2']:
		count_col = 2    # Read count column in Kraken
		level_col = 3    # Level column in Kraken
		tax_col = 4      # Taxonomy column in Kraken
	elif args.filetype == 'krona':
		count_col = 0    # Count column in Krona
		tax_col = 1      # Taxonomy column in Krona
	
	# Read the file
	for line in f:
		try:
			parts = line.strip().split('\t')
			if len(parts) <= max(count_col, tax_col, level_col if level_col != -1 else 0):
				continue  # Skip if not enough columns
			
			# Skip lines with non-numeric count values or comments
			if not parts[count_col].strip().isdigit() or parts[0].startswith('#'):
				continue
				
			count = float(parts[count_col])
			
			# Only process entries matching the requested taxonomy level
			if args.level != 'all' and level_col != -1:
				if not parts[level_col].startswith(args.level):
					continue
			
			if count > 0:
				n.append(count)
				
		except ValueError as e:
			sys.stderr.write(f"Warning: Could not convert count '{parts[count_col] if len(parts) > count_col else 'unknown'}' to float - skipping line\n")
			continue
		except Exception as e:
			sys.stderr.write(f"Warning: Error processing line: {line.strip()} - {str(e)}\n")
			continue
	
	f.close()
	
	# Ensure we have at least some valid data
	if len(n) == 0:
		print("No valid abundance data found for taxonomy level", args.level)
		sys.exit(1)
	
	# calculations
	N = sum(n)  # total number of individuals
	S = len(n)  # total number of species
	
	# calculate all the pi's
	p = []  # store pi's 
	D = 0
	for i in n:  # go through each species
		if i != 0:  # there should not be any zeros
			p.append(i/N)  # pi is the ni/N
			D += i*(i-1)
	
	D = D/(N*(N-1)) if N > 1 else 0
	
	# find the indicated alpha
	if args.value == 'Sh':  # calculate shannon's diversity
		shannons_alpha(p)
	elif args.value == 'BP':  # calculate berger-parker's dominance index
		berger_parkers_alpha(p)
	elif args.value == 'Si':  # calculate Simpson's alpha 
		simpsons_alpha(D)
	elif args.value == 'ISi':  # calculate Inverse Simpson's alpha 
		inverse_simpsons_alpha(D)
	elif args.value == 'F':  # calculate fisher's alpha
		print("Fisher's alpha...loading")
		global N_f
		N_f = sum(n)
		global S_f
		S_f = len(n)
		fishers_alpha()
	else:
		print("Not a supported alpha diversity metric")

if __name__ == "__main__":
    main()

