from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import subprocess

#define function to deliver alignment fragments of a given record based on positions
# corresponding to the n fragments that are gonna be created
def mysplit(a, n):
    k, m = divmod(len(a), n)
    tmp = []
    for i in list(range(n)):
        tmp.append(a[i*k+min(i, m):(i+1)*k+min(i+1, m)])
    return(tmp)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--aln', help="alignment to be split", required=True)
    parser.add_argument('--fnr',help="Number of fragments to split your alignment into",required=True)
    parser.add_argument('--dir', help= "output directory to write alignments to", required=True)
    parser.set_defaults()
    args = parser.parse_args()

    print(args)

    #create the output directory
    subprocess.call(["mkdir", args.dir])
    #loop for opening one file per alignment fragment
    myfiles = []
    for f in list(range(int(args.fnr))):
        tmp = f + 1
        #format number in filename to be able to sort the files by name
        if int(args.fnr) > 0 and int(args.fnr) < 100:
            tmp = f"{tmp:02}" #this format file number as 01 if we have less than 100 fragments
        elif int(args.fnr) >= 100 and int(args.fnr) < 1000:
            tmp = f"{tmp:03}" #This format file number as 001 if we have from 100 to 999 fragments
        file = open(args.dir + "/run_alignment_no_resis.{}.fas".format(tmp), 'a')
        myfiles.append(file)  
    
    j=1 #variable used for creating the list of fragments start positions just once
    #All this block takes the sequences in the alignment one by one, cuts it in n parts, and write
    #the parts to their correspondent files.
    with open (args.aln,"r") as f:
        for record in SeqIO.parse(f, "fasta"):
            i=0 #we reset this variable for every new record we process
            seg = mysplit(record.seq,int(args.fnr)) #keep all alignment fragments in seg
            if j==1:
                #create a file with the starting position of every file to be able to 
                # add it later on when calculating the real position of a mutation in R
                # script transform_split_utationannotation
                with open("mysplits.txt","w") as ms:
                    k, m = divmod(len(record), int(args.fnr)) #divmod is also used in mysplit function
                    mysplits = [0] #0 is the starting position of alignment in file 1
                    for x in list(range(int(args.fnr))):
                        mysplits.append((x+1)*k+min(x+1, m)) #calculates the postion for each of the n files
                    ms.write(str(mysplits))
                    j=0 #change j to 0 to avoid repeating this operation
            #Now we write the fragment of the present record to its correspondent file
            for file in myfiles:
                newrecord = SeqRecord(seg[i],
                               id = record.id,
                                description="")
                SeqIO.write(newrecord, file, "fasta")
                i = i+1 #i controls track of the part of the alignment to be written
    
    for f in myfiles: #we close all files
        f.close()
    
if __name__ == "__main__":
    main()