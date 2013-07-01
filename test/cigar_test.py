# test cigar strings

PROGNAM = "../src/smalt"
FNAM_REF = "cigar_ref.fa.gz"
FNAM_READ1 = "cigar_read1.fq"
FNAM_READ2 = "cigar_read2.fq"

TMPFIL_PREFIX = "TMPcig"

KMER = 13
NSKIP = 2

def smalt_index(df,index_name, fasta_name, kmer, nskip):
    from sys import exit
    from subprocess import call

    tup = (PROGNAM, 'index',
           '-k', '%i' % (int(kmer)),
           '-s', '%i' % (int(nskip)),
           index_name,
           fasta_name)
    df.call(tup, "when indexing")

def smalt_map(df, oufilnam, indexnam, readfil, matefil, typ="fastq", flags=[]):
    from sys import exit
    from subprocess import call
 
    tup = [PROGNAM, 'map']
    if len(flags) > 0:
        tup.extend(flags)
    tup.extend([
           '-f', typ,
           '-o', oufilnam,
           indexnam,
           readfil, matefil])
    df.call(tup, "when mapping")

if __name__ == '__main__':
    from testdata import DataFiles
    
    df = DataFiles()
    
    refnam = df.joinData(FNAM_REF)
    readnamA = df.joinData(FNAM_READ1)
    readnamB = df.joinData(FNAM_READ2)
    indexnam = df.addIndex(TMPFIL_PREFIX)
    oufilnam = df.addTMP(TMPFIL_PREFIX + ".sam")
    
    smalt_index(df,indexnam, refnam, KMER, NSKIP)
    smalt_map(df,oufilnam, indexnam, readnamA, readnamB, "sam", ["-x"])
    
    #print "Test ok."
    
    df.cleanup()
    exit()
