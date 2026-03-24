import os
import pandas as pd
import argparse
import re



df= pd.read_csv(r"C:\Eli\VNTR\muc2exon\69\MUC2_69_done.csv")

outputname = "MUC2_69_done_frequency"

#stack all of the columns
dftall = df.stack()

#turn it back into a dataframe, as it was a panda series once stacked
dftall = dftall.to_frame()

#Rename the one column to 'Stacked'
dftall.columns =['Stacked']

#Find the frequency of each motif in all of the genomes
ranking = dftall['Stacked'].value_counts()

frequencyname = outputname + ("_frequency.csv")
ranking.to_csv(frequencyname, header=False)