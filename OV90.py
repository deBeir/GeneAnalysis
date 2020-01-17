# imports
from Bio import Entrez, SeqIO, AlignIO, ExPASy, Medline, Entrez, SwissProt
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from openpyxl import Workbook

Entrez.email = "email@example.com"

class Search:
    def __init__(self,name, term, keywords):
        self.name = name
        self.term = term
        self.keywords = str(keywords)

    def articlesSearch(self):
        ''' Class method that searches and returns articles that 
            match certain given keywords. The keyword terms must
            be of type 'str'.'''

        search = Entrez.esearch(db="pubmed", term=self.keywords, usehistory="y") 
        record_search = Entrez.read(search)
        print(record_search['Count'])
        fetch_handle = Entrez.efetch(db='pubmed', rettype = 'medline', webenv = record_search["WebEnv"],
                                      Query_Key = record_search['QueryKey'])

        filename = open(self.name+'_data.txt', 'w')
        data = fetch_handle.read()
        filename.write(data) #stores all information inside a file.

        with open("filename.txt","w") as handle:
            records = Medline.parse(handle) #parses all records contained inside the file 
            for record in records: #prints out the title and abstract if the record 
                print('TITLE:', record['TI'], '\n')
                print('ABSTRACT:', record['AB'])

class RetrieveInfo:

    def __init__(self, name, term, start, stop):
        self.name = name
        self.term = str(term)
        self.start = start
        self.stop = stop

    def retrieveGenes(self):
        '''Retrieves the desired sequence and stores it into a
           file (genbank format).'''

        handle = Entrez.efetch(db="nucleotide", id=self.term, rettype="gb", \
                               seq_start=self.start, seq_stop=self.stop)
        record = SeqIO.read(handle, "gb")
        handle.close()

        try:
            self.filename = record.id + ".gb"  #creates file to store the gene data
            try:
                with open(self.filename, 'w') as output_handle: #writes the data into the file created above
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
        '''Retrieves the most important features and stores
           them inside a 'xlsx' file.'''

        wb = Workbook()
        ws = wb.active

        self.record = SeqIO.read(self.filename, "genbank")

        CDS_data = []  #contains all the features collected from all of the CDSs
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

                data.append(gene_id[7:])  #only the iD number
                data.append(self.term)  #NC number
                data.append(product)  #product obtained
                data.append(protein)  #protein iD
                data.append(translation)  #protein sequence

                CDS_data.append(data)

        header = ["GENE_ID", "NCBI_ID", "PRODUCT", "PROTEIN", "TRANSLATION"]
        ws.append(header) #header of the excel file.

        try:
            for data in CDS_data:
                ws.append(data)
        except IOError:
            print("Could not add the information to the file.")

        wb.save(self.name+"_info.xlsx")  #saves the features in an excel file.

        
class Blast:
    def __init__(self, name, term, start, stop):
        self.name = name
        self.term = str(term)
        self.start = start
        self.stop = stop

    def seqBlast(self,file):
        '''Executes a 'blastn' program to search nucleotide databases using
           a nucleotide query. Results are stored inside an 'xml' file.'''

        record = SeqIO.read(file,format = "fasta") 
        result_handle = NCBIWWW.qblast("blastn","nt",record.seq, \
                                      hitlist_size = 10, megablast = True)

        blast_result = open(self.name+"_blast.xml","w")
        blast_result.write(result_handle.read())
        blast_result.close()
        result_handle.close()

    def parseBlast(self,file):
        '''Parses all the sequences obtained in seqBlast and returns useful
           information in case their e-value is less than or equal to 0.05.
           One could set the e-value to 0 if they wish to find the most sig-
           nificant results. '''

        f = open(file)
        blast_records = NCBIXML.parse(f) #parses all sequences contained inside f.

        e_value = 0.05 #or 0, depending on what you want your threshold to be.

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect <= e_value: #prints data only if the e-value is less than or equal to the established threshold
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
        '''Reads single alignment file and prints out.'''

        alignment = AlignIO.read(file, format = str(format)) #single alignment
        print(alignment)

    def multipleAL(self, file, format):
        '''Multiple sequence alignment input/output as alignment objects. '''
        
        for alignment in AlignIO.parse(file,format, seq_count=2): # reads a file with at least 2 alignments
            print ("Alignment length %i" % alignment.get_alignment_length())
            for record in alignment:
                print("%s-%s"%(record.seq, record.id)) 
            print("")

class Clustal:
    def __init__(self, name, term):
        self.name = name
        self.term = str(term)

    def clustalW(self, file):
        '''Will generate alignment given a file as input.
           The file format must be 'fasta' for it to work.'''

        clustal_result = ClustalwCommandline("clustalw2", infile=file)
        print(clustal_result)

        
class Phylo:
    def __init__(self, name, term):
        self.name = name
        self.term = str(term)

    def phylo(self, file, format):
        '''Draws the file phylogram. '''
        tree = Phylo.read(file, format = format) #parses and returns one tree from the given file.
        print(tree)

        Phylo.draw_ascii(tree) #prints an ascii-art rooted phylogram
        tree.rooted = True
        Phylo.draw(tree) # displays a rooted phylogram [rooted = True]

class SearchP:
    def __init__(self,name,term):
        self.name = name
        self.term = str(term)

    def ExPASy(self, accNumber ):
        ''' Given the protein accession number, gets relevant info
            about that protein, like location, function,etc.''' 

        handle = ExPASy.get_sprot_raw(accNumber) #gets a text handle to a raw SwissProt entry at ExPASy.
        record = handle.read()
        print(record)

        
class ProteinBlast:
    def __init__(self, name, term):
        self.name = name
        self.term = str(term)

    def seqBlast(self,file):
        '''Executes a 'blastp' program to search protein databases using
           a protein query. Results are stored inside an 'xml' file in or-
           der to be parsed.'''

        record = SeqIO.read(file,format = "fasta")
        result_handle = NCBIWWW.qblast("blastp","nt",record.seq, \
                                      hitlist_size = 10, megablast = True)

        blast_result = open(self.name+"_blastp.xml","w") #creates file to store the blastp results.
        blast_result.write(result_handle.read())
        blast_result.close()
        result_handle.close()

    def parseBlast(self,file):
        '''Parses all the sequences obtained in seqBlast and returns useful
           information in case their e-value is less than or equal to 0.05.
           One could set the e-value to 0 if they wish to find the most sig-
           nificant results.'''

        f = open(file)
        blast_records = NCBIXML.parse(f)

        e_value = 0.05 # e-value must be at most 0.05 for the results to be significant. 
                       # The closer the e-value is to zero, the more significant.
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
