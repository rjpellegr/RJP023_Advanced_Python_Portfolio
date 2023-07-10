# RJP023_Advanced_Python_Portfolio
A repo of code I learned as part of Louisiana Tech's advanced computational biology (BISC 450C 384) coursework. For use as both a reference and as a submission for my final project for summer quarter 2023.

## Sequence Objects

A basic overview of biopython's Seq objects.

```python
#import necessary libraries/functions
from Bio.Seq import Seq
```


```python
#create a sequence object - format is similar to creating a normal string
my_seq = Seq("GATCG")
```


```python
#sequence objects share many other similarities with strings - in fact, they can easily be converted to a string using the str function
str(my_seq)
```




    'GATCG'




```python
#they can be used in combination with the string interpolation operator %
#for example, to provide comments or labels before the sequence as part of the FASTA format
fasta_format_string = ">Name\n%s" % my_seq
print(fasta_format_string)
```

    >Name
    GATCG



```python
#sequence objects can also be concatenated
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
print("Sequence 1: %s\nSequence 2: %s" % (seq1, seq2))
print("Sequence 1 plus Sequence 2:", seq1 + seq2)
print("Sequence 2 plus Sequence 1:", seq2 + seq1)
```

    Sequence 1: ACGT
    Sequence 2: AACCGG
    Sequence 1 plus Sequence 2: ACGTAACCGG
    Sequence 2 plus Sequence 1: AACCGGACGT



```python
#like a normal string, each character in a sequence object is indexed
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
#we can manipulate these indices using many of the same methods as strings, such as checking its length
print(len(my_seq))
```

    5



```python
#printing values at a specific index
print(my_seq[0])
```

    G



```python
#or counting elements of a specific value
Seq("AAAA").count("AA")
```




    2




```python
#let's use a longer sequence object for further demonstration - this one has 32 nucleotides
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
#as with a string, sequence objects can be sliced by specifying a start and a stop
#to return just the first half of the sequence, we'd use
my_seq[0:15]
```




    Seq('GATCGATGGGCCTAT')




```python
#an optional third parameter, a stride, can also be introduced to skip over certain values
#the syntax is [start:stop:stride] - all values are inclusive!
#to return every third nucleotide starting with index 0 and ending at index 31, we'd use
my_seq[0:31:3]
```




    Seq('GCTGTAGTAAG')




```python
#the default starting index is 0, the default stop index is the last index in the sequence, and the default stride is 1
#the code above can thus be simplified as
my_seq[::3]
```




    Seq('GCTGTAGTAAG')




```python
#using a negative value will start the count from the last index and work backwards
my_seq[::-1]
```




    Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG')




```python
#the join method can be used to join all items in an iterable into a single string, using another string - or a sequence object, in our case - as a separator
#the syntax for this is x.join(y), where x is the separator and y is the iterable
iterable = ("Item1","Item2","Item3")
"separator".join(iterable)
```




    'Item1separatorItem2separatorItem3'




```python
#one potential application for this is assembling contigs into a single sequence
#some sequencers will use "N" to denote gaps in the scaffold or to indicate that a base call could not be made at a given position due to sequence ambiguity, so that will be our spacer
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
spacer = Seq("N" * 10)
#the final assembled sequence will thus look something like this
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
#we can use the upper or lower methods to force all characters in a sequence into the same case
dna_seq = Seq("acgtACGT")
print(dna_seq)
print(dna_seq.upper())
print(dna_seq.lower())
```

    acgtACGT
    ACGTACGT
    acgtacgt



```python
#this is useful when case sensitivity is an issue
"GTAC" in dna_seq
```




    False




```python
"GTAC" in dna_seq.upper()
```




    True




```python
#these string methods can be used for many practical purposes
#for example, we can calculate the length and GC content of a primer - let's use the same 32 bp sequence as before
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
gcCount = 100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
print("Primer: %s" % my_seq)
print("This primer is", len(my_seq), "base pairs long and contains", gcCount, "percent GC content.")
```

    Primer: GATCGATGGGCCTATATAGGATCGAAAATCGC
    This primer is 32 base pairs long and contains 46.875 percent GC content.



```python
#however, Biopython has many built-in functions to make things even easier
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
gcCount = 100 * gc_fraction(my_seq)
print("Primer: %s" % my_seq)
print("This primer is", len(my_seq), "base pairs long and contains", gcCount, "percent GC content.")
```

    Primer: GATCGATGGGCCTATATAGGATCGAAAATCGC
    This primer is 32 base pairs long and contains 46.875 percent GC content.



```python
#other useful Biopython features include the complement method
print("5'- %s - 3'" % my_seq)
print("3'- %s - 5'" % my_seq.complement())
```

    5'- GATCGATGGGCCTATATAGGATCGAAAATCGC - 3'
    3'- CTAGCTACCCGGATATATCCTAGCTTTTAGCG - 5'



```python
#and reverse complement
print("5'- %s - 3'" % my_seq)
print("5'- %s - 3'" % my_seq.reverse_complement())
```

    5'- GATCGATGGGCCTATATAGGATCGAAAATCGC - 3'
    5'- GCGATTTTCGATCCTATATAGGCCCATCGATC - 3'



```python
#The sequence object has a built-in method that transcribes DNA into RNA
#let's start with a coding strand of DNA
coding_dna = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
```


```python
#this transcribe method simply replaces all instances of thymine with uracil
messenger_rna = coding_dna.transcribe()
print("coding DNA:   ", coding_dna)
print("messenger RNA:", messenger_rna)
```

    coding DNA:    ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
    messenger RNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG



```python
#keep in mind that actual biological transcription is a two-step process - it starts with a template strand of DNA and transcribes the mRNA based on its reverse complement
template_dna = Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')
coding_dna = template_dna.reverse_complement()
messenger_rna = coding_dna.transcribe()
print("template DNA: ", template_dna)
print("coding DNA:   ", coding_dna)
print("messenger RNA:", messenger_rna)
```

    template DNA:  CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
    coding DNA:    ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
    messenger RNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG



```python
#the back_transcribe method can be used to convert mRNA back into a piece of coding DNA, which can then be used to find the original template DNA
messenger_rna = Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
coding_dna = messenger_rna.back_transcribe()
template_dna = coding_dna.reverse_complement()
print("messenger RNA:", messenger_rna)
print("coding DNA:   ", coding_dna)
print("template DNA: ", template_dna)
```

    messenger RNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
    coding DNA:    ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
    template DNA:  CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT



```python
#Biopython can translate RNA into into the corresponding amino acid sequence
#stop codons are denoted by an asterisk
messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
amino_acids = messenger_rna.translate()
print("messenger RNA:      ", messenger_rna)
print("amino acid sequence:", amino_acids)
```

    messenger RNA:       AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
    amino acid sequence: MAIVMGR*KGAR*



```python
#to make things even simpler, it can also translate a coding strand of DNA directly
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
amino_acids = coding_dna.translate()
print("coding DNA:         ", coding_dna)
print("amino acid sequence:", amino_acids)
```

    coding DNA:          ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
    amino acid sequence: MAIVMGR*KGAR*



```python
#as can be seen above, this sequence contains multiple stop codons
#you can use the to_stop parameter to stop translation at the first stop codon - this parameter is exclusive, so the stop symbol is not included in the final translation
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
amino_acids = coding_dna.translate(to_stop = True)
print("coding DNA:         ", coding_dna)
print("amino acid sequence:", amino_acids)
```

    coding DNA:          ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
    amino acid sequence: MAIVMGR



```python
#or the stop_symbol parameter to change the asterisk to something else
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
amino_acids = coding_dna.translate(stop_symbol = "!")
print("coding DNA:         ", coding_dna)
print("amino acid sequence:", amino_acids)
```

    coding DNA:          ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
    amino acid sequence: MAIVMGR!KGAR!



```python
#non-standard genetic codes - for example, mitochondrial sequences - may translate codons differently
#in cases like this, you can use the table parameter to specify which codon table is in use
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
amino_acids = coding_dna.translate(table = "Vertebrate Mitochondrial")
print("coding DNA:         ", coding_dna)
print("amino acid sequence:", amino_acids)
```

    coding DNA:          ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
    amino acid sequence: MAIVMGRWKGAR*



```python
#you can also use the NCBI table number instead
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
amino_acids = coding_dna.translate(table = 2)
print("coding DNA:         ", coding_dna)
print("amino acid sequence:", amino_acids)
```

    coding DNA:          ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
    amino acid sequence: MAIVMGRWKGAR*



```python
#the standard start codon is AUG and encodes for methionine
#given a complete coding sequence (CDS), the default translation method will thus usually result in a methionine at the beginning of the amino acid sequence
#however, sometimes you have a complete CDS with a non-standard start codon - a common occurence in bacteria, for example
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
print(gene.translate(table="Bacterial", to_stop = True))
```

    VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDHGWWKQHYEWRGNRWHLHGPPPPPRHHKKAPHDHHGGHGPGKHHR



```python
#the first codon, GTG, has been translated as valine
#while GTG does normally encode for valine, it is also a valid start codon in the bacterial genome and should thus be translated as methionine in this case
#to fix this, we need to specify that the given sequence is a complete CDS using the cds parameter
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
print(gene.translate(table="Bacterial", to_stop = True, cds = True))
```

    MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDHGWWKQHYEWRGNRWHLHGPPPPPRHHKKAPHDHHGGHGPGKHHR



```python
#speaking of tables, the NCBI codon tables used internally for the sequence object's translation method can be visualized with print if needed
from Bio.Data import CodonTable
standard_table = CodonTable.unambiguous_dna_by_id[1]
mito_table = CodonTable.unambiguous_dna_by_id[2]
```


```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
#there are also some other useful properties of these tables such as searching for start or stop codons
print("start codons:", mito_table.start_codons)
print("stop codons: ", mito_table.stop_codons)
```

    start codons: ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
    stop codons:  ['TAA', 'TAG', 'AGA', 'AGG']



```python
#and finding the amino acid encoded by a given codon
print("codon ACG encodes for:", mito_table.forward_table["ACG"])
```

    codon ACG encodes for: T



```python
#like a string, sequence objects are immutable, meaning that their contents can't be edited
#if you want to alter the contents of a sequence, you can use a mutable sequence instead
from Bio.Seq import MutableSeq
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
mutable_seq = MutableSeq(my_seq)
print(mutable_seq)
```

    GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA



```python
#you can also create a mutable sequence object directly from a string in the same way you'd create a normal sequence object
mutable_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
print(mutable_seq)
```

    GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA



```python
#mutable sequence objects can be manipulated in ways that sequence objects cannot
mutable_seq[5] = "C"
print(mutable_seq)
mutable_seq.remove("T")
print(mutable_seq)
mutable_seq.reverse()
print(mutable_seq)
```

    GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA
    GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA
    AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG



```python
#and if you want to make it immutable again, simply make a new sequence object
new_seq = Seq(mutable_seq)
print(new_seq)
```

    AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG



```python
#finally, if for some reason you don't want to work with sequence objects, some Bio.Seq functions also accept normal strings
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
print(reverse_complement(my_string))
print(transcribe(my_string))
print(back_transcribe(my_string))
print(translate(my_string))
```

    CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC
    GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG
    GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG
    AVMGRWKGGRAAG*

## Sequence Annotation

An introduction to Biopython's sequence annotation objects.

```python
#import necessary libraries/functions
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
```


```python
#a SeqRecord object allows you to annotate a Seq object for later use
simple_seq = Seq("GATC")
simple_seq_r = SeqRecord(simple_seq)
```


```python
#if you didn't include an id, name, or description, they will be set as strings by default and can be edited after the fact
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
simple_seq_r.id = "AC12345"
simple_seq_r.name = "testseq"
simple_seq_r.description = "Short test sequence"
print(simple_seq_r)
```

    ID: AC12345
    Name: testseq
    Description: Short test sequence
    Number of features: 0
    Seq('GATC')



```python
#if needed, the annotations dictionary attribute can be used to add other information that doesn't fit under the other categories
simple_seq_r.annotations["evidence"] = "https://youtu.be/iOVbAmknKUk"
print(simple_seq_r)
```

    ID: AC12345
    Name: testseq
    Description: Short test sequence
    Number of features: 0
    /evidence=https://youtu.be/iOVbAmknKUk
    Seq('GATC')



```python
#you can also assign annotations per letter in a sequence, which is very useful when assigning sequencing quality scores for example
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
print(simple_seq_r.seq)
print(simple_seq_r.letter_annotations["phred_quality"])
```

    GATC
    [40, 40, 38, 30]



```python
#let's create a record for a FASTA file containing a complete sequence for Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
from Bio import SeqIO
record = SeqIO.read("NC_005816.fna", "fasta")
print(record)
```

    ID: gi|45478711|ref|NC_005816.1|
    Name: gi|45478711|ref|NC_005816.1|
    Description: gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence
    Number of features: 0
    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')



```python
#the SeqIO function automatically reads and interprets the FASTA format into a sequence record object
#we can look at each of its attributes individually
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
#let's use the seqIO function to make a record from a GenBank file of the same sequence
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
from Bio import SeqIO
record = SeqIO.read("NC_005816.gb", "genbank")
print(record)
```

    ID: NC_005816.1
    Name: NC_005816
    Description: Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence
    Database cross-references: Project:58037
    Number of features: 41
    /molecule_type=DNA
    /topology=circular
    /data_file_division=BCT
    /date=21-JUL-2008
    /accessions=['NC_005816']
    /sequence_version=1
    /gi=45478711
    /keywords=['']
    /source=Yersinia pestis biovar Microtus str. 91001
    /organism=Yersinia pestis biovar Microtus str. 91001
    /taxonomy=['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacteriales', 'Enterobacteriaceae', 'Yersinia']
    /references=[Reference(title='Genetics of metabolic variations between Yersinia pestis biovars and the proposal of a new biovar, microtus', ...), Reference(title='Complete genome sequence of Yersinia pestis strain 91001, an isolate avirulent to humans', ...), Reference(title='Direct Submission', ...), Reference(title='Direct Submission', ...)]
    /comment=PROVISIONAL REFSEQ: This record has not yet been subject to final
    NCBI review. The reference sequence was derived from AE017046.
    COMPLETENESS: full length.
    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')



```python
#any PROJECT or DBLINK links are stored in the dbxrefs list
record.dbxrefs
```




    ['Project:58037']




```python
#entries in the features table are automatically stored as SeqFeature objects
len(record.features)
```




    41




```python
#SeqFeature objects are useful for providing information about a specified region of a sequence
#for example, we can create a location with ambiguous end points
from Bio import SeqFeature
#we start by defining the start and end points
start_pos = SeqFeature.AfterPosition(5)
end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
#and then create the location proper with the SimpleLocation method
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
print(my_location)
```

    [>5:(8^9)]



```python
#we can look at each location individually
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
#since they're subclasses of integers, we can also use them as integers
int(my_location.start)
```




    5




```python
int(my_location.end)
```




    9




```python
#for an exact location, you can simply pass numbers to the constructors
exact_location = SeqFeature.SimpleLocation(5,9)
print(exact_location)
```

    [5:9]



```python
#the parameters are stored as ExactPosition objects
exact_location.start
```




    ExactPosition(5)




```python
#SeqRecord objects can be sliced and added together
#let's start with our genbank file
record = SeqIO.read("NC_005816.gb", "genbank")
print(record)
```

    ID: NC_005816.1
    Name: NC_005816
    Description: Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence
    Database cross-references: Project:58037
    Number of features: 41
    /molecule_type=DNA
    /topology=circular
    /data_file_division=BCT
    /date=21-JUL-2008
    /accessions=['NC_005816']
    /sequence_version=1
    /gi=45478711
    /keywords=['']
    /source=Yersinia pestis biovar Microtus str. 91001
    /organism=Yersinia pestis biovar Microtus str. 91001
    /taxonomy=['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacteriales', 'Enterobacteriaceae', 'Yersinia']
    /references=[Reference(title='Genetics of metabolic variations between Yersinia pestis biovars and the proposal of a new biovar, microtus', ...), Reference(title='Complete genome sequence of Yersinia pestis strain 91001, an isolate avirulent to humans', ...), Reference(title='Direct Submission', ...), Reference(title='Direct Submission', ...)]
    /comment=PROVISIONAL REFSEQ: This record has not yet been subject to final
    NCBI review. The reference sequence was derived from AE017046.
    COMPLETENESS: full length.
    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')



```python
#this is a circular genome, so the origin of the sequence isn't so clear-cut
#if we wanted to shift the origin for whatever reason, we can slice the SeqRecord at two places and add them together like so
shifted = record[2000:] + record[:2000]
```


```python
#the length stays the same
print(len(record))
print(len(shifted))
```

    9609
    9609



```python
#but the origin has now been shifted 2000 base pairs over
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
#however, slicing is selective in what annotations are preserved
#in this example, we've lost a feature
len(record.features)
```




    41




```python
len(shifted.features)
```




    40




```python
#specifically, the database cross references
record.dbxrefs
```




    ['Project:58037']




```python
shifted.dbxrefs
```




    []




```python
#to preserve this, we need to explicitly assign the annotations in question
shifted.dbxrefs = record.dbxrefs[:]
shifted.annotations = record.annotations.copy()
```


```python
#the dbxrefs have now been restored
shifted.dbxrefs
```




    ['Project:58037']




```python
#the reverse_complement method from Seq object also works on SeqRecords
record = SeqIO.read("NC_005816.gb", "genbank")
rc = record.reverse_complement(id="TESTING")
```


```python
#however, it should be noted that the id, name, description, annotations, and dbxrefs are not transferred by default
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    NC_005816.1 9609 41 1 13
    TESTING 9609 41 0 0

