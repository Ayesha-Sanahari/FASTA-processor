"""
Author : Ayesha Sanahari
Date: 16/02/2021
Python script to process sequences in FASTA format (FASTA processor)
"""

class fasta:
    # Constructor method with input parameters
    def __init__(self, file_name):
        self.file_name = file_name

    # A method to remove unwanted/unknown characters from a FASTA nucleotide sequence
    def remove_unwanted_characters_from_NA(file_name):
        from Bio import SeqIO
        import re

        handle = open(file_name, "r")
        sequences = {}

        for record in SeqIO.parse(handle, "fasta"):
            # replace any non-GATC characters in your sequences with nothing
            sequence = re.sub('[^GATC]', "", str(record.seq).upper())
            sequence = sequence.replace("\n", "")
            sequence = sequence.replace("\r", "")
            sequence = sequence.replace("N", "")
            sequences[sequence] = record.description

        output_file = open("remove_unwanted_characters_from_NA.fasta", "w+")
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n\n")
        output_file.close()
        handle.close()


    # A method to remove unwanted/unknown characters from a FASTA amino acid sequence
    def remove_unwanted_characters_from_AA(file_name):
        from Bio import SeqIO
        import re

        handle = open(file_name, "r")
        sequences = {}

        for record in SeqIO.parse(handle, "fasta"):
            # replace any non-GATC characters in your sequences with nothing
            sequence = re.sub('[^KNTRSIMQHPLEDAGVOYCWF]', "", str(record.seq))
            sequence = sequence.replace("\n", "")
            sequence = sequence.replace("\r", "")
            sequences[sequence] = record.description

        output_file = open("remove_unwanted_characters_from_AA.fasta", "w+")
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n\n")
        output_file.close()
        handle.close()


# ************


class SingleFasta(fasta):
    def __init__(self, file_name):
        # Using the super method to pass object variables to the parent class
        super().__init__(file_name)

    #• A method to calculate and append the length of an input FASTA sequence to the respective FASTA header.
    def addLengthToHeader(file_name):
        from Bio import SeqIO

        record = SeqIO.read(file_name, "fasta")
        # length = len(record.seq)
        # len =str(length)

        # write the fasta seq with length in header into a fasta file
        with open('FASTA_file_with_length.fasta', 'w') as output:
            output.write('>' + record.description + ' length:' + str(len(record.seq)) + '\n')
            output.write(str(record.seq))

    # *************

    # •	A method to calculate and append the AT or GC content of an input FASTA sequence and append it to the FASTA header.
    # The FASTA sequence, sequence type, and content type (AT or GC) should be the input parameters for the method.
    def get_ATorGC_content(file_name, seq_type, content_type):
        from Bio import SeqIO

        record = SeqIO.read(file_name, "fasta")
        content = ''
        sequence = record.seq
        length = len(record.seq)

        if seq_type == "DNA":
            if content_type == "AT":
                A_count = sequence.upper().count('A')
                T_count = sequence.upper().count('T')
                AT_content = (A_count + T_count) / length

                data_AT = '>' + record.description + " AT_content: " + str(
                    round(AT_content, 3)) + '\n' + sequence + '\n\n'
                content += data_AT
                return content

                # write the fasta seq with AT or GC Content in header into a fasta file
                with open('FASTA_file_with_ATorGC_content.fasta', 'w') as output:
                    output.write('>' + record.description + " AT_content: " + str(round(AT_content, 3)) + '\n' + str(sequence) + '\n\n')

        if seq_type == "RNA" or "mRNA":
            if content_type == "AT":
                A_count = sequence.upper().count('A')
                T_count = sequence.upper().count('U')
                AT_content = (A_count + T_count) / length

                data_AU = '>' + record.description + " AT_content: " + str(
                    round(AT_content, 3)) + '\n' + sequence + '\n\n'
                content += data_AU
                return content

                # write the fasta seq with AT or GC Content in header into a fasta file
                with open('FASTA_file_with_ATorGC_content.fasta', 'w') as output:
                    output.write('>' + record.description + " AT_content: " + str(round(AT_content, 3)) + '\n' + str(
                        sequence) + '\n\n')

        if seq_type == "DNA" or "RNA" or "mRNA":
            if content_type == "GC":
                G_count = sequence.upper().count('G')
                C_count = sequence.upper().count('C')
                GC_content = (G_count + C_count) / length

                data_GC = '>' + record.description + " GC_content: " + str(
                    round(GC_content, 3)) + '\n' + sequence + '\n\n'
                content += data_GC
                return content

                # write the fasta seq with AT or GC Content in header into a fasta file
                with open('FASTA_file_with_ATorGC_content.fasta', 'w') as output:
                    output.write('>' + record.description + " GC_content: " + str(round(GC_content, 3)) + '\n' + str(sequence) + '\n\n')


# ************

class MultiFasta(fasta):
    def __init__(self, file_name):
        # Using the super method to pass object variables to the parent class
        super().__init__(file_name)
        self.Seq_count = self.count_no_of_fasta_Seqs(file_name)

    # A method to get seq count in a multiple fasta file
    def count_no_of_fasta_Seqs(self, file_name):

        from Bio.SeqIO.FastaIO import SimpleFastaParser

        count = 0
        total_len = 0
        with open(file_name) as in_handle:
            for title, seq in SimpleFastaParser(in_handle):
                count += 1
                total_len += len(seq)

        return ("%i records with total sequence length %i" % (count, total_len))

    # ***************

    # A method to split multiple FASTA sequences in a single file to several FASTA files.
    def splitMultipleFasta(file_name):
        from Bio import SeqIO
        record_iterator = SeqIO.parse(open(file_name), "fasta")

        # for i, batch in enumerate(batch_iterator(record_iter, 1000)):
        for i, batch in enumerate(record_iterator, 0):
            filename = "Sequence_%i.fasta" % (i + 1)
            with open(filename, "w") as handle:
                count = SeqIO.write(batch, handle, "fasta")
            print("Wrote %i records to %s" % (count, filename))

    # *****************

    #• A method to calculate and append the length of each sequence to the respective FASTA header for a given input file containing multiple FASTA sequences.
    def addLengthToHeader(file_name):
        from Bio import SeqIO
        content = ''
        for seq_record in SeqIO.parse(file_name, "fasta"):
            length = len(seq_record)
            seq = seq_record.seq
            aaa = '>' + seq_record.description + ' length: ' + str(length) + '\n' + seq + '\n\n'
            content += aaa
        print(content)

        # write the fasta sequence with length in header into a fasta file
        with open('multiple_fasta_with_length.fasta', 'w') as output:
            output.write(str(content))

    # ***********

    # A method to combine individual FASTA files in a certain folder into a single FASTA file.
    def combine_many_FASTA_files(DIR):
        import os

        # instead of DIR, you can enter the directory/path (a folder containing that all fasta files you want to combine)
        # DIR = 'fasta_files/'
        with open("single_fasta_file.fasta", 'w') as output:
            for file in os.listdir(DIR):
                file_holder = open(os.path.join(DIR, file))
                for line in file_holder:
                    output.write(line)
                file_holder.close()


# *************
class Sequence:
    # Constructor method with input parameters
    def __init__(self, sequence):
        self.sequence = sequence


class DNAseq(Sequence):
    def __init__(self, sequence):
        # Using the super method to pass object variables to the parent class
        super().__init__(sequence)
        self.AT_content = self.get_ATcontent(sequence)
        self.GC_content = self.get_GCcontent(sequence)

    #•	Methods to calculate AT and GC content in DNA sequences
    def get_ATcontent(self, sequence):
        length = len(sequence)
        A_count = sequence.upper().count('A')
        T_count = sequence.upper().count('T')
        AT_content = (A_count + T_count) / length
        return round(AT_content, 3)

    def get_GCcontent(self, sequence):
        length = len(sequence)
        G_count = sequence.upper().count('G')
        C_count = sequence.upper().count('C')
        GC_content = (G_count + C_count) / length
        return round(GC_content, 3)


class RNAseq(Sequence):
    def __init__(self, sequence):
        # Using the super method to pass object variables to the parent class
        super().__init__(sequence)
        self.AT_content = self.get_ATcontent(sequence)
        self.GC_content = self.get_GCcontent(sequence)

    # •	Methods to calculate AT and GC content in RNA sequences
    def get_ATcontent(self, sequence):
        length = len(sequence)
        A_count = sequence.upper().count('A')
        T_count = sequence.upper().count('U')
        AT_content = (A_count + T_count) / length
        return round(AT_content, 3)

    def get_GCcontent(self, sequence):
        length = len(sequence)
        G_count = sequence.upper().count('G')
        C_count = sequence.upper().count('C')
        GC_content = (G_count + C_count) / length
        return round(GC_content, 3)



# Implementation
if __name__ == '__main__':
    # creating the first object
    file_obj1 = MultiFasta("multiple_fasta_file.fasta")
    file_obj2 = SingleFasta("group_12.fasta")
    seq_obj1 = DNAseq("AAATCGCTAGCTTAGGGACCAT")
    seq_obj2 = RNAseq("AUUGCAUCCAUUGUACCAAUAUGCCU")

    file_obj5 = MultiFasta("Together all cacao proteins_.fasta")
    file_obj6 = MultiFasta("Together all Ara proteins_.fasta")


print("Seq count in multi fasta file:", file_obj1.Seq_count)
print(MultiFasta.splitMultipleFasta("multiple_fasta_file.fasta"))
print(MultiFasta.addLengthToHeader("multiple_fasta_file.fasta"))
print(MultiFasta.combine_many_FASTA_files('fasta_files/'))
print(SingleFasta.addLengthToHeader("group_12.fasta"))
print("AT_content:", seq_obj1.AT_content)
print("AT_content:", seq_obj2.AT_content)
print("GC_content:", seq_obj1.GC_content)
print("GC_content:", seq_obj2.GC_content)

print(SingleFasta.get_ATorGC_content("group_11.fasta", "DNA", "AT"))
print(SingleFasta.get_ATorGC_content("group_11.fasta", "DNA", "GC"))

print(fasta.remove_unwanted_characters_from_NA("multiple_fasta_file_with_unwanted_letters.fasta"))
print(fasta.remove_unwanted_characters_from_AA("OSDREB_aa_sequences_with_unwanted_letters.fasta"))


print("Theobroma cacao")
print("Seq count in multi fasta Together all cacao proteins file:", file_obj5.Seq_count)
print("Seq count in multi fasta Together all Ara proteins file:", file_obj6.Seq_count)
