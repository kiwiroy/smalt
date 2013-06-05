# test diff string printing

PROGNAM = "./diffstr_test"

DIFFFSTR_DATA = [("13S65M3I4M2X12M3X20M12S", "13S65M3I41M12S", "M 65 I 3 M 41"),
                 ("20M3D77M4X8M", "20M3D89M", "M 20 D 3 M 89"),
                 ("2S124M16X46M2X8M1S", "2S196M1S", "M 196"),
                 ("1S61M2X1M1I62M1X2M", "1S64M1I65M", "M 64 I 1 M 65"),
                 ("61M", "61M", "M 61")
                 ]
import re

CIGSTR = re.compile("\nCIGAR_SSAHA\s+([MID\s\d]+)\s*\nCIGAR_EXT\s+(\S+)\s*\nCIGAR_SAM\s+(\S+)\s*\nCIGAR_SAMX\s+(\S+)\s*$")
#CIGSTR = re.compile("^CIGAR_SAM", re.MULTILINE)
def checkDiffStr(diffstr):
    from sys import exit
    from subprocess import check_output, CalledProcessError 

    
    tup = (PROGNAM, diffstr[0])

    try:
        rc = check_output(tup)
    except CalledProcessError as dfserr:
        exit("%s %s returned %i" % (PROGNAM, diffstr[0], dfserr.returncode))


    m = CIGSTR.search(rc)
    if m:
        if m.group(4) != diffstr[0]:
            exit("DIFF strings do not agree: '%s' vs '%s'" % (m.group(4), diffstr[0]))
        if m.group(3) != diffstr[1]:
            exit("DIFF strings do not agree: '%s' vs '%s'" % (m.group(3), diffstr[1]))        
        if m.group(2) != diffstr[1]:
            exit("DIFF strings do not agree: '%s' vs '%s'" % (m.group(2), diffstr[1]))
        if m.group(1).strip() != diffstr[2]:
            exit("DIFF strings do not agree: '%s' vs '%s'" % (m.group(1), diffstr[2]))         
    else:
        print "Not matched\n**%s**\n" % rc

    
if __name__ == '__main__':

    for dat in DIFFFSTR_DATA:
        checkDiffStr(dat)
    
    exit(0)
    
