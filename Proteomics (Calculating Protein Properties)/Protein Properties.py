# Author: Rajiv Snape
# Project Title: Protein Properties Calculation

# These are the packages and modules needed to carry my analysis and culation
# SeqIO to read in fasta files
# ProteinAnalysis to perform biochemical stat calculatiosn
# arr to use arrys
# matplot to plot
# csv to handle csv files

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import array as arr
import matplotlib.pyplot as plt
import csv

# Declaring an array for the isoelectric values
Iso_array = arr.array('f')

# Declaration of header string
header = "Seq ID, Mol Weight, Iso Point, Gravy \n"
# This for loops parses the proteome file with the help of the SeqIO module
# each record is a protein within the file
for record in SeqIO.parse("Influenza Proteome.fasta", "fasta"):
    print(record.id)
    MyProteins = ProteinAnalysis(str(record.seq))
    
    Mol_weight = MyProteins.molecular_weight()
    print(Mol_weight)

    Iso_points = MyProteins.isoelectric_point()
    print(Iso_points)

    Gravy_score = MyProteins.gravy()
    print(Gravy_score)
    print("\n")

    # adds Isolelectric values to the declared array
    Iso_array.append(Iso_points)

    # data variable careated to store and record the calculations
    data = f"{record.id}, {Mol_weight}, {Iso_points}, {Gravy_score} \n"
    # compling all the data together
    header = header + data
print(header)

# creating csv file
f = open("Protein Calculation.csv", 'w')
# writing to the csv file
f.write(header)
f.close()

# For loop to print values in the array
for i in (Iso_array):
    print (i, end =" ")
print()

# use matplot to create as histogram of the values on the created array.
plt.hist(Iso_array, bins=None, range=None, density=True, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, data=None)
plt.show()











