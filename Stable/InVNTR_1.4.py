from Bio.Seq import Seq
import os
import pandas as pd
import argparse
import re

#create a variable for any errors
error = "Errors:\n"

#argument parser for commandline use.
parser = argparse.ArgumentParser(description='Extract VNTR')
parser.add_argument('-f', '--folder', type=str, metavar='', required=True, help='Path to folder. ex: E:\HPRC')
parser.add_argument('-l', '--length', type=int, metavar='', nargs='?', const=1, default= '0', help='length of consensus motif. ex: 30 : DO NOT USE IF USING DELIMITER, LENGTH TAKES PRIORITY')
parser.add_argument('-d', '--delimiter', type=str, metavar='', nargs='?', const=1, default= '0', help='consistent beginning of motif: DO NOT USE IF USING LENGTH, LENGTH TAKES PRIORITY')
parser.add_argument('-s', '--start', type=str, metavar='', required=True, help='beginning of VNTR. ex: GCTA')
parser.add_argument('-e', '--end', type=str, metavar='', nargs='?', const=1, default= 'n', help='end of VNTR. ex: GCTA')
parser.add_argument('-ne', '--no_end', type=int, metavar='', nargs='?', const=1, default= '10000', help='If there is no end, when should the program cut out? default is 10000')
parser.add_argument('-mal', '--max_allele_length', type=int, metavar='', nargs='?', const=1, default= '10000', help='Maximum nucleotide length that will be written to the allele file output, per file. Default is 10000')
parser.add_argument('-b', '--before', type=str, metavar='', nargs='?', const=1, default= 'n', help='Is the start value of the VNTR before the VNTR? y OR n')
parser.add_argument('-a', '--after', type=str, metavar='', nargs='?', const=1, default= 'n', help='Is the end value of the VNTR after the VNTR? y OR n')
parser.add_argument('-t', '--type', type=str, metavar='', nargs='?', const=1, default= 0, help='Filetype ex: .fasta (default is .fa or fasta)')
parser.add_argument('-n', '--name', type=str, metavar='', nargs='?', const=1, default= 'VNTR', help='output name)')
args = parser.parse_args()

#Lets define the user inputted variables first:


#Where are the files?
folder = args.folder

#What is the length of the motif?
motif_length = int(args.length)

#is there a delimiter?
if motif_length == 0:

#What is the delimiter?
    delimiter = str(0)
    delimiter = str(args.delimiter)
    
    if delimiter == '0':
        print("ERROR: neither length nor delimiter specified.")
        error += ("neither length nor delimiter specified \n")


#What is the beginning of the VNTR, and is it the beginning of the VNTR or directly BEFORE the VNTR?
beginning = args.start #consistent beginning of the VNTR
before = args.before # n means it is NOT before the VNTR (it is within the VNTR). y means it is before the VNTR.

#What is the end of the VNTR, and is it the end of the VNTR or directly AFTER the VNTR?
ending = args.end #consistent ending of the VNTR
after = args.after # n means it is not AFTER the VNTR (it is within the VNTR). y means it is after the VNTR.

#What filetype?

filetype = args.type
if filetype == 0:
    filetype = (".fa", ".fasta")

#What name should be appended to the output file?
outputname = args.name

#make a variable the length of the beginning
beg_length = len(beginning)

#If the before variable is NOT n, then add the length of the beginning to the index of the beginning to not include what is before the VNTR
if before.find('n') != -1:
    beg_adjustment = 0
    endrc_adjustment = beg_length
else:
    beg_adjustment = beg_length
    endrc_adjustment = 0

#measure the length of the ending variable: this will help us change the index value for the location of the ending if it is within the VNTR
end_length = len(ending)

if end_length == 1:
    ne = args.no_end

#If the after variable is n, then add the length of the ending to the index of the end to encapsulate the entire VNTR
if after.find('n') != -1:
    end_adjustment = end_length
    begrc_adjustment = 0
else:
    end_adjustment = 0
    begrc_adjustment = end_length

#set the reverse complement of the beginning and ending units
seqbeg = Seq(beginning)
seqend = Seq(ending)

#set the reverse complement of the first two repeat units as the end
beginningrc = str(seqend.reverse_complement()) #reverse complement of end to beginning
endingrc = str(seqbeg.reverse_complement()) #reverse complement of beginning to end

#create the dataframe, which will ultimately be written as a csv
df = pd.DataFrame()

#set the index from 0 - 10000, 10000 will be the maximum number of repeat units that can be output
dfindex = [*range(1, 10000, 1)]
df[0] = pd.Series(dfindex)
df = df.set_index(0)

#create a variable for the the final index size, this will be the number that the index is cut down to at the end of the program
finalindexsize = 2

#create max allele length variable (this variable controls how long of an allele to save to the allele file per file)
max_allele_length = args.max_allele_length

#create allele variable
alleles = ""

def find_config_end(string):
    for index, char in enumerate(string):
        if char in ['C', 'T', 'G', 'A'] and char.isupper():
            return index
    return -1  # Return -1 if no match is found

#create allele length variable which will create a csv which with the allele length of each file
allele_length = ""

if end_length == 1: #if the end length variable was not specified, and therefore this has NO END sequence, and will at default run for 10000 bp
    for filename in os.listdir(folder): #for loop which opens each genome to a variable, 'x', and the filename to a variable, 'filename'
        if filename.endswith(filetype):
            with open(os.path.join(folder, filename)) as x:
                print(filename)

                #reads the genome from variable 'x' to variable 'a'
                a = x.read()

                #remove line breaks (\n) and save to variable 'b'
                b = a.replace("\n","")
                
                #search for and make an index with the location of the beginning, if it can be found
                beg_index = b.find(beginning)

                #was the beginning of the repeat found? If not beg_index will = -1, otherwise it will be a positive number
                if beg_index < 0: #in this case, the sequence is likely in the reverse complement

                    #find the locations of the beginning of the repeat in the reverse complement
                    end_indexrc = b.find(endingrc)

                    if end_indexrc < 0: #cannot find beginning in forward or reverse
                        error += filename
                        error += (" has no sign of tandem repeat in forward, reverse \n")
                        
                    else: 
                        #adjust for whether or not the start and end points are before or after the VNTR
                        endrc = end_indexrc+endrc_adjustment

                        #lets save the reverse complement of the tandem repeat to new variable 'c'
                        c = b[endrc-ne:endrc]

                        #lets see if there is a contig break somewhere in this
                        c_index = c.rfind('>')

                        if c_index > 0: #if there is a contig break we're going to only keep everything in the first contig
                            error += filename #add filename to error string
                            error += (" had multiple contigs, kept only the first \n") #add error message
                            configname100 = c[c_index:c_index+100] #from the > marking the beginning of the config name, lets take the next 100 characters
                            configendindex = find_config_end(configname100) #find the first capital ATCG after the >
                            c = c[c_index+configendindex:] #lets replace c with only that first capital ATCG

                        #lets make this reverse complement a Seq object 'd'

                        d = Seq(c)

                        #now we can reverse complement it and save it as a string

                        #add the filename to the allele variable
                        alleles += filename
                        alleles += "\n"

                        tandemrepeat = str(d.reverse_complement())
                        alleles += tandemrepeat[:max_allele_length] # keep only the max allele length
                        alleles += "\n\n"

                        #add the filename to the allele length variable
                        allele_length += filename
                        allele_length += ","
                        #find the length of the allele length and save it to variable allele length
                        allele_length += str(len(tandemrepeat[:max_allele_length]))
                        allele_length += "\n"


                        if motif_length == 0:
                            #split the tandemrepeat variable into motifs by the delimiter
                            rx =  [delimiter+e for e in tandemrepeat.rsplit(delimiter) if e]

                        if motif_length != 0:
                            #split tandemrepeat variable into motifs by motif_length from left to right.
                            rx = [tandemrepeat[idx:idx + motif_length] for idx in range(0, len(tandemrepeat), motif_length)]

                        #append to beginning of list so that it doesn't overwrite the first motif when I turn it into a pd series
                        rx.insert(0, "overwrite me")

                        #update the final index size by comparing the current length of repeat unit to the current length of the finalindexsize
                        if len(rx) > finalindexsize:
                            finalindexsize = len(rx) 
                            
                        #add the filename as the header of a new column with the repeats below
                        df[filename] = pd.Series(rx)
                        
                else:
                    beg = beg_index+beg_adjustment

                    #put it into a variable
                    tandemrepeat = b[beg:beg+ne]

                    tandemrepeat_index = tandemrepeat.find('>') #lets see if there is a contig break somewhere in this

                    if tandemrepeat_index > 0: #if there is a contig break we're going to only keep everything in the first contig 
                        error += filename
                        error += (" had multiple contigs, kept only the first \n")
                        tandemrepeat = tandemrepeat[:tandemrepeat_index] #keep only from beginning to >, the first contig


                    #add the filename to the allele variable
                    alleles += filename
                    alleles += "\n"

                    #add the allele to the alleles variable
                    alleles += tandemrepeat[:max_allele_length] # keep only the max allele length
                    alleles += "\n\n"     

                    #add the filename to the allele length variable
                    allele_length += filename
                    allele_length += ","

                    #find the length of the allele length and save it to variable allele length
                    allele_length += str(len(tandemrepeat[:max_allele_length]))
                    allele_length += "\n"

                    #split the tandemrepeat variable into motifs by the delimiter

                    if motif_length == 0:
                        #split the tandemrepeat variable into motifs by the delimiter
                        rx =  [delimiter+e for e in tandemrepeat.rsplit(delimiter) if e]


                    if motif_length != 0:
                        #split tandemrepeat variable into motifs by motif_length from left to right.
                        rx = [tandemrepeat[idx:idx + motif_length] for idx in range(0, len(tandemrepeat), motif_length)]

                    #append to beginning of list so that it doesn't overwrite the first motif when I turn it into a pd series
                    rx.insert(0, "overwrite me")

                    #update the final index size by comparing the current length of repeat unit to the current length of the finalindexsize
                    if len(rx) > finalindexsize:
                        finalindexsize = len(rx) 
                        
                    #add the filename as the header of a new column with the repeats below
                    df[filename] = pd.Series(rx)

else: #giant forloop for searching with an end
    for filename in os.listdir(folder): #for loop which opens each genome to a variable, 'x', and the filename to a variable, 'filename'
        if filename.endswith(filetype):
            with open(os.path.join(folder, filename)) as x:
                print(filename)

                #reads the genome from variable 'x' to variable 'a'
                a = x.read()

                #remove line breaks (\n) and save to variable 'b'
                b = a.replace("\n","")

                #search for and make an index with the location of the beginning, if it can be found
                beg_index = b.find(beginning)

                #was the beginning of the repeat found? If not beg_index will = -1, otherwise it will be a positive number
                if beg_index < 0: #in this case, the sequence is likely in the reverse complement

                    #find the locations of the beginning of the repeat in the reverse complement
                    beg_indexrc = b.find(beginningrc)

                    if beg_indexrc < 0: #cannot find beginning in either forward or reverse... can we find the ends?
                        
                        #make an index with the location of the forward ending
                        end_index = b.find(ending)

                        if end_index < 0: #cannot find ending in forward

                            #find the locations of the end of the repeat in the reverse complement
                            end_indexrc = b.find(endingrc)

                            if end_indexrc < 0: #cannot find any sign of the tandem repeat
                                error += filename
                                error += (" has no sign of tandem repeat in forward, reverse \n")

                            else: #cannot find reverse beginning
                                error += filename
                                error += (" has only reverse end \n")
                                
                        else: #cannot find forward beginning
                            error += filename
                            error += (" has only forward end \n")

                    else:
                        #find the locations of the end of the repeat in the reverse complement
                        end_indexrc = b.find(endingrc)

                        if end_indexrc < 0: #cannot find reverse end
                            error += filename
                            error += (" has only reverse beginning \n")
                        
                        else: 
                            #adjust for whether or not the start and end points are before or after the VNTR
                            begrc = beg_indexrc+begrc_adjustment
                            endrc = end_indexrc+endrc_adjustment

                            if endrc < begrc:
                                error += filename
                                error += (" Your provided start was found after your provided end, maybe switch them? \n")
                                print("Your provided start was found after your provided end, maybe switch them?")

                            #lets save the reverse complement of the tandem repeat to new variable 'c'
                            c = b[begrc:endrc]


                            #lets see if there is a contig break somewhere in this
                            c_index = c.rfind('>')

                            if c_index > 0: #if there is a contig break we're going to only keep everything in the first contig
                                error += filename #add filename to error string
                                error += (" had multiple contigs, kept only the first \n") #add error message
                                configname100 = c[c_index:c_index+100] #from the > marking the beginning of the config name, lets take the next 100 characters
                                configendindex = find_config_end(configname100) #find the first capital ATCG after the >
                                c = c[c_index+configendindex:] #lets replace c with only that first capital ATCG

                            #lets make this complement a Seq object 'd'
                            d = Seq(c)


                            #add the filename to the allele variable
                            alleles += filename
                            alleles += "\n"

                            #now we can complement it and save it as a string
                            tandemrepeat = str(d.reverse_complement())
                            alleles += tandemrepeat[:max_allele_length]
                            alleles += "\n\n"

                            #add the filename to the allele length variable
                            allele_length += filename
                            allele_length += ","
                            #find the length of the allele length and save it to variable allele length
                            allele_length += str(len(tandemrepeat[:max_allele_length]))
                            allele_length += "\n"

                            if motif_length == 0:
                                #split the tandemrepeat variable into motifs by the delimiter
                                rx =  [delimiter+e for e in tandemrepeat.rsplit(delimiter) if e]

                            if motif_length != 0:
                                #split tandemrepeat variable into motifs by motif_length from left to right.
                                rx = [tandemrepeat[idx:idx + motif_length] for idx in range(0, len(tandemrepeat), motif_length)]

                            #append to beginning of list so that it doesn't overwrite the first motif when I turn it into a pd series
                            rx.insert(0, "overwrite me")

                            #update the final index size by comparing the current length of repeat unit to the current length of the finalindexsize
                            if len(rx) > finalindexsize:
                                finalindexsize = len(rx) 
                                
                            #add the filename as the header of a new column with the repeats below
                            df[filename] = pd.Series(rx)
                            

                else:
                    #make an index with the location of the ending
                    end_index = b.find(ending)

                    if end_index < 0: #cannot find forward end
                        error += filename
                        error += (" has only forward beginning \n")

                    else:
                        #adjust for whether or not the start and end points are before or after the VNTR
                        beg = beg_index+beg_adjustment
                        end = end_index+end_adjustment

                        if end < beg:
                                error += filename
                                error += (" Your provided start was found after your provided end, maybe switch them? \n")
                                print("Your provided start was found after your provided end, maybe switch them?")


                        #put it into a variable
                        tandemrepeat = b[beg:end]

                        tandemrepeat_index = tandemrepeat.find('>') #lets see if there is a contig break somewhere in this

                        if tandemrepeat_index > 0: #if there is a contig break we're going to only keep everything in the first contig
                            error += filename
                            error += (" had multiple contigs, kept only the first \n")
                            tandemrepeat = tandemrepeat[:tandemrepeat_index] #keep only from beginning to >, the first contig


                        #add the filename to the allele variable
                        alleles += filename
                        alleles += "\n"

                        #add the allele to the alleles variable
                        alleles += tandemrepeat[:max_allele_length]
                        alleles += "\n\n"     

                        #add the filename to the allele length variable
                        allele_length += filename
                        allele_length += ","

                        #find the length of the allele length and save it to variable allele length
                        allele_length += str(len(tandemrepeat[:max_allele_length]))
                        allele_length += "\n"

                        #split the tandemrepeat variable into motifs by the delimiter

                        if motif_length == 0:
                            #split the tandemrepeat variable into motifs by the delimiter
                            rx =  [delimiter+e for e in tandemrepeat.rsplit(delimiter) if e]


                        if motif_length != 0:
                            #split tandemrepeat variable into motifs by motif_length from left to right.
                            rx = [tandemrepeat[idx:idx + motif_length] for idx in range(0, len(tandemrepeat), motif_length)]

                        #append to beginning of list so that it doesn't overwrite the first motif when I turn it into a pd series
                        rx.insert(0, "overwrite me")

                        #update the final index size by comparing the current length of repeat unit to the current length of the finalindexsize
                        if len(rx) > finalindexsize:
                            finalindexsize = len(rx) 
                            
                        #add the filename as the header of a new column with the repeats below
                        df[filename] = pd.Series(rx)

#drop all rows below the finalindexsize, ensuring that the size of the dataframe is confined to where we actually have repeats
df = df.drop(labels=range(finalindexsize, 10000), axis=0)

csvname = outputname + (".csv")
#save the dataframe to a csv 'out.csv'
df.to_csv(csvname)

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

#Add a Frequency category of the text file
#error += ("\nFrequency:\n")

#Add the ranking to the frequency section, as a string
#error += ranking.to_string()

errorsname = outputname + ("_errors.txt") 
#save errors and frequency variable as textfile
with open(errorsname, "w") as text_file:
    text_file.write(error)

allelename = outputname + ("_alleles.txt")
#save alleles variable as textfile
with open(allelename, "w") as text_file:
    text_file.write(alleles)

allele_length_name = outputname + ("_allele_length.csv")
with open(allele_length_name, "w") as text_file:
    text_file.write(allele_length)

#To do... figure out how to deal with VNTR's with motifs of variable length.
