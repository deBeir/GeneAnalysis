import OV90

# gene STX4 used as an example
# missing tests for classes Align, Clustal and Phylo


def test_RetrieveInfo():

    STX4 = OV90.RetrieveInfo('STX4', 'NC_000016.10',31033095,31040168)

    STX4.retrieveGenes()
    STX4.retrieveAnnotations()
    STX4.retrieveFeatures()

def test_Blastn():

    STX4 = OV90.Blast('STX4','NC_000016.10',31033095,31040168)
    STX4.seqBlast()

def test_Align():
    STX4 = OV90.Align('STX4','NC_000016.10')

    STX4.singleAL()
    STX4.multipleAL()

def test_SearchP():

    STX4 = OV90.SearchP('STX4','NC_000016.10')
    STX4.ExPASy('Q12846')

def test_Blastp():
    STX4 = OV90.ProteinBlast('STX4', 'NC_000016.10')
    STX4.seqBlast('STX4.fasta')
    STX4.seqParse('STX4.xml')

if __name__ == "__main__":

    #test_RetrieveInfo()
    #test_Blastn()
    #test_Align()
    #test_SearchP()
    #test_Blastp()





