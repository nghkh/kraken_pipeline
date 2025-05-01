#!/usr/bin/env python
################################################################
# beta_diversity.py - Calculates Bray-Curtis dissimilarity metrics
# between multiple community samples
################################################################
import os, sys, argparse
import operator
from time import gmtime
from time import strftime
import numpy as np

####################################################################
# Main method
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Calculate Bray-Curtis dissimilarity between multiple samples.")
    parser.add_argument('-i','--input','--input-files',
        '--inputs', required=True, dest='in_files', nargs='+',
        help='Input files (one per community) for which to compare for \
        bray-curtis dissimiliarity metrics')
    parser.add_argument('--type', required=False, default='single', dest='filetype',
        choices=['single','simple','bracken','kreport','kreport2','krona'],
        help='Type of input file[s]: single, simple [tab-delimited, specify --cols], \
            bracken, kreport, kreport2, krona. See docs for details')
    parser.add_argument('--cols','--columns', dest='cols', required=False, default='1,2',
        help='Specify category/counts separated by single comma: cat,counts (1 = first col)')
    parser.add_argument('--level', '-l', dest='lvl', required=False, default='all',
        choices=['all', 'D', 'P', 'C', 'O', 'F', 'G', 'S'],
        help='Taxonomy level for which to compare samples. Default: all')
    args = parser.parse_args()

    #################################################
    # Test input files
    in2counts = {}
    if args.filetype == 'single' and len(args.in_files) > 1:
        sys.stderr.write("Please specify only one file for '--type single'\n")
        exit(1)
    for f in args.in_files:
        if not os.path.isfile(f):
            sys.stderr.write("File %s not found\n" % f)
            exit(1)

    #################################################
    # Determine columns for extracting
    categ_col = -1
    count_col = -1
    if args.filetype in ['single','simple']:
        if ',' not in args.cols:
            sys.stderr.write("Please specify column as 'a,b' where a = column of category, \
            b = column of first count\n")
            exit(1)
        else:
            [categ_col, count_col] = args.cols.split(',')
            if not categ_col.isdigit():
                sys.stderr.write("%s is not an integer\n" % categ_col)
                exit(1)
            elif not count_col.isdigit():
                sys.stderr.write("%s is not an integer\n" % count_col)
                exit(1)
            categ_col = int(categ_col) - 1
            count_col = int(count_col) - 1
    elif args.filetype == "bracken":
        categ_col = 0  # Taxonomy column
        count_col = 5  # Abundance column
        taxlvl_col = 2  # Level column
    elif args.filetype in ["kreport", "kreport2"]:
        categ_col = 4  # Name column
        count_col = 2  # Reads column
        taxlvl_col = 3  # Level column
    elif args.filetype == "krona":
        categ_col = 1  # Taxonomy column
        count_col = 0  # Count column
    
    #################################################
    # STEP 1: READ IN SAMPLES
    i2totals = {}       # Total counts per sample
    i2counts = {}       # Counts per category per sample
    i2names = {}        # Sample names
    num_samples = 0     # Number of samples
    num_categories = 0  # Number of categories
    
    if args.filetype == "single":
        # ALL SAMPLE COUNTS WITHIN A SINGLE FILE
        header = True
        i_file = open(args.in_files[0],'r')
        for line in i_file:
            l_vals = line.strip().split("\t")
            # Read header
            if header:
                s_count = 0
                for i in range(count_col, len(l_vals)):
                    i2names[s_count] = l_vals[i]
                    i2counts[s_count] = {}
                    i2totals[s_count] = 0
                    s_count += 1
                num_samples = s_count
                header = False
            else:
                # Otherwise, save counts
                s_count = 0
                curr_categ = l_vals[categ_col]
                for i in range(count_col, len(l_vals)):
                    if int(l_vals[i]) > 0:
                        i2totals[s_count] += int(l_vals[i])
                        i2counts[s_count][curr_categ] = int(l_vals[i])
                    s_count += 1
                num_categories += 1
        i_file.close()
    else:  # For bracken, kraken, kraken2 and krona
        num_samples = 0
        i2names = {}
        i2totals = {}
        i2counts = {}
        taxonomies = {}  # Store taxonomy information

        for f in args.in_files:
            i_file = open(f,'r')
            i2names[num_samples] = f
            i2totals[num_samples] = 0
            i2counts[num_samples] = {}
            taxonomies[num_samples] = {}

            for line in i_file:
                l_vals = line.strip().split("\t")

                # Skip empty lines, non-numeric counts, or comments
                if len(l_vals) <= max(count_col, categ_col, taxlvl_col if 'taxlvl_col' in locals() else -1):
                    continue
                if not l_vals[count_col].isdigit() or l_vals[0].startswith('#'):
                    continue

                count = int(l_vals[count_col])
                if count <= 0:
                    continue
                    
                if args.filetype == "krona":
                    # Handle special case for Krona format
                    for i in range(count_col, len(l_vals)):
                        taxonomy_prefix = args.lvl.lower() + "__"
                        if (l_vals[i].startswith(taxonomy_prefix) or args.lvl == "all"):
                            tax_col = i
                            tax_name = l_vals[tax_col]
                            i2totals[num_samples] += count
                            if tax_name not in i2counts[num_samples]:
                                i2counts[num_samples][tax_name] = 0
                            i2counts[num_samples][tax_name] += count
                elif args.filetype == "bracken":
                    # Handle bracken format - improved check for taxonomic level
                    # For bracken files with .G.report, .S.report, etc. in the filename
                    # the level is included in the filename but the internal level might be different
                    if args.lvl == "all" or (taxlvl_col < len(l_vals) and 
                                           (l_vals[taxlvl_col] == args.lvl or 
                                            l_vals[taxlvl_col][0] == args.lvl)):
                        tax_id = l_vals[categ_col]
                        taxonomies[num_samples][tax_id] = l_vals[0]
                        i2totals[num_samples] += count
                        if tax_id not in i2counts[num_samples]:
                            i2counts[num_samples][tax_id] = 0
                        i2counts[num_samples][tax_id] += count
                else:
                    # Handle Kraken report formats
                    if args.lvl == "all" or (taxlvl_col < len(l_vals) and l_vals[taxlvl_col][0] == args.lvl):
                        tax_id = l_vals[categ_col]
                        taxonomies[num_samples][tax_id] = l_vals[0]
                        i2totals[num_samples] += count
                        if tax_id not in i2counts[num_samples]:
                            i2counts[num_samples][tax_id] = 0
                        i2counts[num_samples][tax_id] += count

            i_file.close()
            
            # Check if we found any data for this sample
            if i2totals[num_samples] == 0:
                sys.stderr.write(f"Warning: No data found in {f} for taxonomy level {args.lvl}\n")
                
            num_samples += 1
            
    # Verify that we have data to process
    if num_samples == 0:
        sys.stderr.write("Error: No valid samples found\n")
        exit(1)
        
    if all(total == 0 for total in i2totals.values()):
        sys.stderr.write(f"Error: No data found for taxonomy level {args.lvl} in any sample\n")
        exit(1)
    
    #################################################
    # STEP 2: CALCULATE BRAY-CURTIS DISSIMILARITIES
    
    bc = np.zeros((num_samples, num_samples))
    for i in range(0, num_samples):
        i_tot = i2totals[i]
        for j in range(i+1, num_samples):
            j_tot = i2totals[j]
            C_ij = 0.0
            # Skip calculation if either sample has zero counts
            if i_tot == 0 or j_tot == 0:
                bc_ij = 1.0  # Maximum dissimilarity when data is missing
            else:
                for cat in i2counts[i]:
                    if cat in i2counts[j]:
                        C_ij += min(i2counts[i][cat], i2counts[j][cat])
                # Calculate Bray-Curtis dissimilarity
                bc_ij = 1.0 - ((2.0 * C_ij) / float(i_tot + j_tot))
            
            bc[i][j] = bc_ij
            bc[j][i] = bc_ij

    #################################################
    # STEP 3: PRINT MATRIX OF BRAY-CURTIS DISSIMILARITIES
    
    # Print sample information
    for i in i2names:
        sys.stdout.write("#%i\t%s (%i reads)\n" % (i, i2names[i], i2totals[i]))
        
    # Print headers
    sys.stdout.write("x")
    for i in range(num_samples):
        sys.stdout.write("\t%i" % i)
    sys.stdout.write("\n")
    
    # Print dissimilarity matrix
    for i in range(num_samples):
        sys.stdout.write("%i" % i)
        for j in range(num_samples):
            if i <= j:
                sys.stdout.write("\t%0.3f" % bc[i][j])
            else:
                sys.stdout.write("\tx.xxx")
        sys.stdout.write("\n")

####################################################################
if __name__ == "__main__":
    main()
