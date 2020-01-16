from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
from Bio import ExPASy
from Bio import Medline
from Bio import Entrez

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SwissProt
from Bio.Align.Applications import ClustalwCommandline

from openpyxl import Workbook


Entrez.email = "email@example.com"


class Search:

    def __init__(self,name, term, keywords):
        self.name = name
        self.term = term
        self.keywords = str(keywords)

    def articlesSearch(self):

        ''' Given certain keywords, returns articles that match said terms of search '''

        search= Entrez.esearch(db="pubmed", term=self.keywords, usehistory="y")
        record_search = Entrez.read(search)
        print(record_search['Count'])

        fetch_handle = Entrez.efetch(db='pubmed', rettype = 'medline', webenv = record_search["WebEnv"],
                                      Query_Key = record_search['QueryKey'])

        filename = open(self.name+'_data.txt', 'w')
        data = fetch_handle.read()
        filename.write(data)

        with open("filename.txt","w") as handle:
            records = Medline.parse(handle)
            for record in records:
                print('TITLE:', record['TI'], '\n')
                print('ABSTRACT:', record['AB'])

class RetrieveInfo:

    def __init__(self, name, term, start, stop):
        self.name = name
        self.term = str(term)
        self.start = start
        self.stop = stop

    def retrieveGenes(self):
        ''' Retrieves a file containing the desired sequence. '''

        handle = Entrez.efetch(db="nucleotide", id=self.term, rettype="gb", \
                               seq_start=self.start, seq_stop=self.stop)

        record = SeqIO.read(handle, "gb")
        handle.close()

        try:
            self.filename = record.id + ".gb"  # creates a file to store the gene data

            try:
                # writes the data into the file created above
                with open(self.filename, 'w') as output_handle:
                    SeqIO.write(record, output_handle, "gb")

            except:
                print("Error! Could not open file.")
        except:
            print("Error! Could not fetch gene info.")

    def retrieveAnnotations(self):
        ''' Retrieves the annotations and external references. '''

        self.record = SeqIO.read(self.filename, "genbank")
        print("#=== ANNOTATIONS FROM GENE {0} ===#\n{1}".format(self.name, self.record.annotations))
        print("\n")
        print("#=== EXTERNAL REFERENCES FROM GENE {0} ===#\n{1}".format(self.name, self.record.dbxrefs))

        return self.record

    def retrieveFeatures(self):
        ''' Retrieves the most important features. '''

        wb = Workbook()
        ws = wb.active

        self.record = SeqIO.read(self.filename, "genbank")

        CDS_data = []  # contains all the features collected from all of the CDSs
        data = []  # features of each CDS

        # for all features of type CDS, returns the gene iD, term, product, protein and translation.
        for feat in range(len(self.record.features)):

            if self.record.features[feat].type == "CDS":

                if len(self.record.features[feat].qualifiers["db_xref"]) == 2:
                    gene_id = self.record.features[feat].qualifiers["db_xref"][1]

                else:
                    gene_id = self.record.features[feat].qualifiers["db_xref"][0]

                product = self.record.features[feat].qualifiers["product"][0]
                protein = self.record.features[feat].qualifiers["protein_id"][0]
                translation = self.record.features[feat].qualifiers["translation"][0]

                data.append(gene_id[7:])  # only the iD number
                data.append(self.term)  # NC number
                data.append(product)  # product obtained
                data.append(protein)  # protein iD
                data.append(translation)  # protein sequence

                CDS_data.append(data)

        header = ["GENE_ID", "NCBI_ID", "PRODUCT", "PROTEIN", "TRANSLATION"]
        ws.append(header)

        try:
            for data in CDS_data:
                ws.append(data)

        except IOError:
            print("Error.")

        wb.save("despair.xlsx")  # saves the features inside an excel file.

class Blast:

    def __init__(self, name, term, start, stop):
        self.name = name
        self.term = str(term)
        self.start = start
        self.stop = stop

    def seqBlast(self,file):

        record = SeqIO.read(file,format = "fasta")
        result_handle = NCBIWWW.qblast("blastn","nt",record.seq, \
                                      hitlist_size = 10, megablast = True)

        blast_result = open(self.name+"_blast.xml","w")
        blast_result.write(result_handle.read())
        blast_result.close()
        result_handle.close()

    def parseBlast(self,file):

        f = open(file)
        blast_records = NCBIXML.parse(f)

        e_value = 0.05 # or 0.0

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect <= e_value:
                        print("Sequence",alignment.title)
                        print("Length",alignment.length)
                        print("e-value",hsp.expect)
                        print(hsp.query[0:75] + '...')
                        print(hsp.match[0:75] + '...')
                        print(hsp.sbjct[0:75] + '...')

class Align:

    def __init__(self,name,term):
        self.name = name
        self.term = str(term)

    def singleAL(self, file, format):

        alignment = AlignIO.read(file, format = str(format))  # reads a single alignment file
        print (alignment)

    def multipleAL(self, file, format):

        for alignment in AlignIO.parse(file,format, seq_count=2):    # reads a file with more than one alignmnet

            print ("Alignment length %i" % alignment.get_alignment_length())

            for record in alignment:
                print("%s-%s"%(record.seq, record.id))

            print("")

class Clustal:

    def __init__(self, name, term):
        self.name = name
        self.term = str(term)

    def clustalW(self, file):

        clustal_result = ClustalwCommandline("clustalw2", infile=file)
        print(clustal_result)

class Phylo:

    def __init__(self, name, term):
        self.name = name
        self.term = str(term)

    def phylo(self, file, format):
        tree = Phylo.read(file, format = format)
        print(tree)

        Phylo.draw_ascii(tree)
        tree.rooted = True
        Phylo.draw(tree)

class SearchP:

    def __init__(self,name,term):
        self.name = name
        self.term = str(term)


    def ExPASy(self, accNumber ):

        handle = ExPASy.get_sprot_raw(accNumber)
        record = handle.read()
        print(record)

class ProteinBlast:

    def __init__(self, name, term, start, stop):
        self.name = name
        self.term = str(term)
        self.start = start
        self.stop = stop

    def seqBlast(self,file):

        record = SeqIO.read(file,format = "fasta")
        result_handle = NCBIWWW.qblast("blastp","nt",record.seq, \
                                      hitlist_size = 10, megablast = True)

        blast_result = open(self.name+"_blastp.xml","w")
        blast_result.write(result_handle.read())
        blast_result.close()
        result_handle.close()

    def parseBlast(self,file):

        f = open(file)
        blast_records = NCBIXML.parse(f)

        e_value = 0.05 # or 0.0

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect <= e_value:
                        print("Sequence",alignment.title)
                        print("Length",alignment.length)
                        print("e-value",hsp.expect)
                        print(hsp.query[0:75] + '...')
                        print(hsp.match[0:75] + '...')
                        print(hsp.sbjct[0:75] + '...')

