# Tests for module sequence.c
# $Id$
#
import sys, os, re

diffstr = re.compile("^Files .+ are identical")

TESTPROGRAMS = ["sequenceReadWrite_test", "sequenceCompress_test", "sequenceDeCompress_test"]
TMPFILA = "tmpA.fa"
TMPFILB = "tmpB.fa"
TMPDIFFIL = "tmpdiff.txt"

FILES = ["sequence_test_1.fa", "sequence_test_2.fa",
         "sequence_test_3.fq", "sequence_test_4.fq"]


def testReadWrite(df, cnam, fnam, tnum):
    from subprocess import call
    from testdata import openFile
    is_ok = False
    tmpfilA = df.addTMP(TMPFILA)
    diffilnam = df.addTMP(TMPDIFFIL)
    
    fcomtup = ["./%s" % cnam, fnam, tmpfilA, "%i" % tnum]
    df.call(fcomtup)

    tup = ["diff", "-aiws", fnam, tmpfilA]
    df.call(tup, oufilnam = diffilnam)
    
    fdiff = openFile(diffilnam, 'r')
    while 1:
        lin = fdiff.readline()
        if not lin: break
        m = diffstr.search(lin)
        if m:
            is_ok = True
            break

    fdiff.close()
    return is_ok

def testCompress(df, cnam, fnam, s_start, s_end):
    from subprocess import call
    from testdata import openFile
    is_ok = False
    tmpfilA = df.addTMP(TMPFILA)
    tmpfilB = df.addTMP(TMPFILB)
    diffilnam = df.addTMP(TMPDIFFIL)
    
    tup = ["./%s" % cnam, fnam, tmpfilA, tmpfilB, "%i" % s_start, "%i" % s_end]
    df.call(tup);
    
    tup = ["diff", "-aiws", tmpfilA, tmpfilB]
    df.call(tup, oufilnam=diffilnam)

    fdiff = openFile(diffilnam, 'r')

    while 1:
        lin = fdiff.readline()
        if not lin: break
        m = diffstr.search(lin)
        if m:
            is_ok = True
            break

    fdiff.close()
    return is_ok

def testSetReadWrite(df):

    testcom = TESTPROGRAMS[0]
    is_ok = True
    for fil in df.datafiles:
        for i in range(0,3):
            if not testReadWrite(df, testcom, fil, i):
                is_ok = False
                break
        if not is_ok: break
    if is_ok: print "Sequence test %s: ok" % testcom
    else: print "Sequence test: %s %i failed on '%s'" % (testcom, i, fil)
    return  is_ok

def testSetCompress(df):
    arglist = [(df.datafiles[0], (61, 312), (2, 28), (0, 2007540), (2007400,3), (2007491,9)),
               (df.datafiles[3], (0,32),(1,32), (7,25), (8,25), (9,25), (17,13), (16,14))]
    

    is_ok = True
    for testcom in TESTPROGRAMS[1:2]:
        for a in arglist:
            for b in a[1:]:
                if not testCompress(df, testcom, a[0], b[0], b[1]):
                    is_ok = False
                    break
            if not is_ok: break
        if is_ok: print "Sequence test %s: ok" % testcom
        else:
            print "Sequence test %s: %s %i %i failed on '%s'" % (testcom, b[0], b[1], a[0])
        if not is_ok: break

    return  is_ok

if __name__ == '__main__':
    from testdata import DataFiles
    df = DataFiles(FILES)

    is_ok = testSetReadWrite(df)
    if is_ok:
        is_ok = testSetCompress(df)
    if is_ok:
        df.cleanup()
     
    sys.exit(not is_ok)
   
