# testing that forward/reverse complement alignment strings
# are consistent

from sys import path
path.append('../misc')

PROGNAM = "../src/smalt"
#PROGNAM = "smalt-0.7.0.1" # this should fail tests
WORKDIR = "./tmp"
TMPFIL_PREFIX = "TMP"

KMER = 11
NSKIP = 1

REFSEQ = ">REF\n" \
         "aaaaaaaataaataggtaattattaaaaggtaacatatatatatatatgtgtatatattt\n"\
         "attatttttttttgttttggtttcattctattccatttttatagcttatccatgctatgt\n"\
         "aattctttcatataccttactttttttttttcatcctttatttatattaaataaaattat\n"\
         "tgttatttgaagattacatttttgtcttttatttttttttcattatatttttacatacat\n"\
         "tttatacttttcatatttgtctattagctttttcaaacctttcattattgttttcattct\n"\
         "ttttttaacgaaaactattcatctcaaaaatataagatattttatatgacgaatgccatt\n"\
         "gtattttttgttacgtaaaacctgacttcttcaaggaaaacacatgcgcattttcaccaa\n"\
         "tttttgcctaagcttattataaaaagtatattaaatgtatgacttgttcaggtgagacaa\n"\
         "gtatagatttaaactattttggctggatattgttatataaagtatattaaatatatgacc\n"\
         "tgttcaggtaacacacatataaacacatacatatatatatatatatatatatatatatat\n"

READ = "@READ\n"\
       "gcgcgcaattctttcatataccttacttttttttttttttcatcgtttatttatattaataaaattatccggccggccggccgg"

# Expect the following alignment:
# READ  gcgcgcaattctttcatataccttacttttttttttttttcatcgtttatttatattaa-taaaattatccggccggccggccgg
# REF   ++++++aattctttcatataccttacttttttttttt---catcctttatttatattaaataaaattat++++++++++++++++
# SAM 121 6S31M3I19M1D9M16S

TARGET_REFPOS = 121
TARGET_CIGARSTR = "6S31M3I19M1D9M16S"
 
def smalt_index(df, index_name, fasta_name, kmer, nskip):
    from subprocess import call
    from sys import exit
    tup = (PROGNAM, 'index',
           '-k', '%i' % (int(kmer)),
           '-s', '%i' % (int(nskip)),
           index_name,
           fasta_name)
    df.call(tup, "smalt index failed")
 
    return 


def smalt_map(df,index_name, fasta_name, mapped_name):
    from subprocess import call
    from sys import exit
    tup = (PROGNAM, 'map',
           '-o', mapped_name,
           index_name,
           fasta_name)
    df.call(tup, "smalt map failed")

    return 

def reverseComplement(df, fastq_name_F, fastq_name_R):
    from subprocess import call
    from sys import exit

    tup = ("./sequenceReverseComplement_test",
           fastq_name_F,
           fastq_name_R)
    df.call(tup, "reverse complement failed")

    return

def checkSAM(filnam_sam):
    from SAM import Sam, openFile
    from sys import exit

    sam = Sam()
    infil = openFile(filnam_sam, 'r')
    sam.next(infil)
    infil.close()

    if not sam.ok:
         exit("ERROR: Could not parse file '%s'" % (filnam_sam))

    if sam.pos != TARGET_REFPOS:
        exit("ERROR: wrong reference position %i (target:%i)" % (sam.pos, TARGET_REFPOS))

    if sam.cigar != TARGET_CIGARSTR:
        exit("ERROR: wrong CIGAR string '%s' (target:'%s')" % (sam.cigar, TARGET_CIGARSTR))
    return

def writeAsFile(filnam, datastr):
    from SAM import openFile
    oufil = openFile(filnam, 'w')
    oufil.write(datastr + '\n')
    oufil.close()
           
def process(df):
    from os import path
    
    filnam_refseq = df.addTMP(TMPFIL_PREFIX + '_ref.fa')
    filnam_index = df.addIndex(TMPFIL_PREFIX + '_idx')
    filnam_read =  df.addTMP(TMPFIL_PREFIX + "_read.fq")
    filnam_RCread = df.addTMP(TMPFIL_PREFIX + "_read.RC.fq")
    filnam_mapped = df.addTMP(TMPFIL_PREFIX + '_mapped.sam')
    filnam_RCmapped = df.addTMP(TMPFIL_PREFIX + '_mapped.RC.sam')

    writeAsFile(filnam_refseq, REFSEQ)
    writeAsFile(filnam_read, READ)
    smalt_index(df,filnam_index, filnam_refseq, KMER, NSKIP)
    smalt_map(df,filnam_index, filnam_read, filnam_mapped)
    reverseComplement(df,filnam_read, filnam_RCread)
    smalt_map(df,filnam_index, filnam_RCread, filnam_RCmapped)
    checkSAM(filnam_mapped)
    checkSAM(filnam_RCmapped)
    
    return

if __name__ == '__main__':
    from subprocess import call
    from sys import exit
    from testdata import DataFiles
    
    df = DataFiles()
    process(df)
    df.cleanup()

    #print "Test of consistent F/R alignment was successful."
    
    exit(0)
    
