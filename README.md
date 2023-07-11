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

## Sequence Input Output
A look at biopython's SeqIO function.

```python
#import necessary libraries/functions
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
```


```python
#to read sequence data in a file, we use the function SeqIO.parse() - this function is an iterator that creates SeqRecord objects
#it requires two arguments: a filename and a format
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk
#to parse the above GenBank file ls_orchid.gbk, we'd use
for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    Z78532.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    Z78531.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    Z78530.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    Z78529.1
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    Z78527.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    Z78525.1
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    Z78524.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    Z78523.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    Z78522.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    Z78521.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    Z78520.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    Z78519.1
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    Z78518.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    Z78517.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    Z78516.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    Z78515.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    Z78514.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    Z78513.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    Z78512.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    Z78511.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    Z78510.1
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    Z78509.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    Z78508.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    Z78507.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    Z78506.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    Z78505.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    Z78504.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    Z78503.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    Z78502.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    Z78501.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    Z78500.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    Z78499.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    Z78498.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    Z78497.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    Z78496.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    Z78495.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    Z78494.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    Z78493.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    Z78492.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    Z78491.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    Z78490.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    Z78489.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    Z78488.1
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    Z78487.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78486.1
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    Z78485.1
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    Z78484.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    Z78483.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    Z78482.1
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    Z78481.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    Z78480.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    Z78479.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    Z78478.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    Z78477.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    Z78476.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    Z78475.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    Z78474.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    Z78473.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    Z78472.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    Z78471.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78470.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    Z78469.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    Z78468.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    Z78467.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    Z78466.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    Z78465.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    Z78464.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    Z78463.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    Z78462.1
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    Z78461.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78460.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    Z78459.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    Z78458.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    Z78457.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    Z78456.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78455.1
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    Z78454.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    Z78453.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    Z78452.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    Z78451.1
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    Z78450.1
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    Z78449.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    Z78448.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    Z78447.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    Z78446.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    Z78445.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    Z78444.1
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    Z78443.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    Z78442.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    Z78440.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
#if you only want to extract specific data, you can also use the parse function as part of a list comprehension
#for example, if we only want the identifiers we can use
identifiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank")]
identifiers
```




    ['Z78533.1',
     'Z78532.1',
     'Z78531.1',
     'Z78530.1',
     'Z78529.1',
     'Z78527.1',
     'Z78526.1',
     'Z78525.1',
     'Z78524.1',
     'Z78523.1',
     'Z78522.1',
     'Z78521.1',
     'Z78520.1',
     'Z78519.1',
     'Z78518.1',
     'Z78517.1',
     'Z78516.1',
     'Z78515.1',
     'Z78514.1',
     'Z78513.1',
     'Z78512.1',
     'Z78511.1',
     'Z78510.1',
     'Z78509.1',
     'Z78508.1',
     'Z78507.1',
     'Z78506.1',
     'Z78505.1',
     'Z78504.1',
     'Z78503.1',
     'Z78502.1',
     'Z78501.1',
     'Z78500.1',
     'Z78499.1',
     'Z78498.1',
     'Z78497.1',
     'Z78496.1',
     'Z78495.1',
     'Z78494.1',
     'Z78493.1',
     'Z78492.1',
     'Z78491.1',
     'Z78490.1',
     'Z78489.1',
     'Z78488.1',
     'Z78487.1',
     'Z78486.1',
     'Z78485.1',
     'Z78484.1',
     'Z78483.1',
     'Z78482.1',
     'Z78481.1',
     'Z78480.1',
     'Z78479.1',
     'Z78478.1',
     'Z78477.1',
     'Z78476.1',
     'Z78475.1',
     'Z78474.1',
     'Z78473.1',
     'Z78472.1',
     'Z78471.1',
     'Z78470.1',
     'Z78469.1',
     'Z78468.1',
     'Z78467.1',
     'Z78466.1',
     'Z78465.1',
     'Z78464.1',
     'Z78463.1',
     'Z78462.1',
     'Z78461.1',
     'Z78460.1',
     'Z78459.1',
     'Z78458.1',
     'Z78457.1',
     'Z78456.1',
     'Z78455.1',
     'Z78454.1',
     'Z78453.1',
     'Z78452.1',
     'Z78451.1',
     'Z78450.1',
     'Z78449.1',
     'Z78448.1',
     'Z78447.1',
     'Z78446.1',
     'Z78445.1',
     'Z78444.1',
     'Z78443.1',
     'Z78442.1',
     'Z78441.1',
     'Z78440.1',
     'Z78439.1']




```python
#you can also use the next function to move through entries step by step
record_iterator = SeqIO.parse("ls_orchid.gbk", "genbank")

first_record = next(record_iterator)
print(first_record.id)
print(first_record.description)

second_record = next(record_iterator)
print(second_record.id)
print(second_record.description)
```

    Z78533.1
    C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Z78532.1
    C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
#in cases where we want to access records in any order, we can use the default Python lists
records = list(SeqIO.parse("ls_orchid.gbk", "genbank"))
print("Found %i records" % len(records))

print("The first record")
first_record = records[0]  # remember, Python counts from zero
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))

print("The last record")
last_record = records[-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    Found 94 records
    The first record
    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
#let's go back to our original SeqRecord object ls_orchid.gbk and look at how we can extract data from it
record_iterator = SeqIO.parse("ls_orchid.gbk", "genbank")
first_record = next(record_iterator)
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
#the annotations dictionary can be printed directly, if desired
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
#as well as the keys
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
#and the values
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78533'], 1, '2765658', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium irapeanum', 'Cypripedium irapeanum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
#to extract a list of species from the ls_orchid.gbk file, we can look in the annotations dictionary under 'source'
>>> print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
#or 'organism'
>>> print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
#using this, let's build a list of species
all_species = []
for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
    all_species.append(seq_record.annotations["organism"])
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
#we can do the same thing with a list comprehension, with the same results
all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank")
]
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
#unfortunately, FASTA files aren't annotated in a standardized way like GenBank files
#this means extracting a species list from a FASTA file is a little more complicated because we have to break up and extract information from the description line
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.fasta
#we'll use the above ls_orchid.fasta file for this example
all_species = []
for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
    all_species.append(seq_record.description.split()[1])
print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
#again, we can also do the same thing more easily with a list comprehension
all_species = [
    seq_record.description.split()[1]
    for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta")
]
print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
#another thing we can do is alter the attributes of a SeqRecord directly
record_iterator = SeqIO.parse("ls_orchid.fasta", "fasta")
first_record = next(record_iterator)
first_record.id
```




    'gi|2765658|emb|Z78533.1|CIZ78533'




```python
#let's change the id
first_record.id = "new_id"
first_record.id
```




    'new_id'




```python
#you can change the way the FASTA file is output by modifying both the id and the description attributes
record_iterator = SeqIO.parse("ls_orchid.fasta", "fasta")
first_record = next(record_iterator)
first_record.id = "new_id"
first_record.description = first_record.id + " | " + "mutations induced randomly"
print(first_record.format("fasta")[:200])
```

    >new_id | mutations induced randomly
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGG



```python
#now that we've looked at sequence input (reading files), it's time to look at sequence output (writing files)
#the SewIO.write function takes three arguments: SeqRecord objects, a handle or filename, and a sequence format
#we'll start with the records
rec1 = SeqRecord(
    Seq(
        "MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
        "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
        "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
        "SSAC",
    ),
    id="gi|14150838|gb|AAK54648.1|AF376133_1",
    description="chalcone synthase [Cucumis sativus]",
)

rec2 = SeqRecord(
    Seq(
        "YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ"
        "DMVVVEIPKLGKEAAVKAIKEWGQ",
    ),
    id="gi|13919613|gb|AAK33142.1|",
    description="chalcone synthase [Fragaria vesca subsp. bracteata]",
)

rec3 = SeqRecord(
    Seq(
        "MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC"
        "EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP"
        "KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN"
        "NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV"
        "SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW"
        "IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT"
        "TGEGLEWGVLFGFGPGLTVETVVLHSVAT",
    ),
    id="gi|13925890|gb|AAK49457.1|",
    description="chalcone synthase [Nicotiana tabacum]",
)
my_records = [rec1, rec2, rec3]
```


```python
#now we can take these records, assign it a name and a format, and export it as a file
SeqIO.write(my_records, "my_example.faa", "fasta")
#when it's done, the write function will also return the number of records written to the file
```




    3




```python
#let's check it out
for seq_record in SeqIO.parse("my_example.faa", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    gi|14150838|gb|AAK54648.1|AF376133_1
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVG...SAC')
    184
    gi|13919613|gb|AAK33142.1|
    Seq('YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAP...WGQ')
    84
    gi|13925890|gb|AAK49457.1|
    Seq('MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKE...VAT')
    389



```python
#we can use the write function in combination with a SeqRecord iterator to convert between sequence file formats
records = SeqIO.parse("ls_orchid.gbk", "genbank")
count = SeqIO.write(records, "my_example.fasta", "fasta")
print("Converted %i records" % count)
```

    Converted 94 records



```python
#to make things even simpler, Biopython also has a built-in convert function
count = SeqIO.convert("ls_orchid.gbk", "genbank", "my_example.fasta", "fasta")
print("Converted %i records" % count)
#keep in mind that if the output file already exists, IT WILL BE OVERWRITTEN WITHOUT WARNING!!!
```

    Converted 94 records



```python
#let's try taking a file of nucleotide sequences and outputting a file containing their reverse complements
for record in SeqIO.parse("ls_orchid.gbk", "genbank"):
    print(record.id)
    print(record.id)
```

    Z78533.1
    Z78533.1
    Z78532.1
    Z78532.1
    Z78531.1
    Z78531.1
    Z78530.1
    Z78530.1
    Z78529.1
    Z78529.1
    Z78527.1
    Z78527.1
    Z78526.1
    Z78526.1
    Z78525.1
    Z78525.1
    Z78524.1
    Z78524.1
    Z78523.1
    Z78523.1
    Z78522.1
    Z78522.1
    Z78521.1
    Z78521.1
    Z78520.1
    Z78520.1
    Z78519.1
    Z78519.1
    Z78518.1
    Z78518.1
    Z78517.1
    Z78517.1
    Z78516.1
    Z78516.1
    Z78515.1
    Z78515.1
    Z78514.1
    Z78514.1
    Z78513.1
    Z78513.1
    Z78512.1
    Z78512.1
    Z78511.1
    Z78511.1
    Z78510.1
    Z78510.1
    Z78509.1
    Z78509.1
    Z78508.1
    Z78508.1
    Z78507.1
    Z78507.1
    Z78506.1
    Z78506.1
    Z78505.1
    Z78505.1
    Z78504.1
    Z78504.1
    Z78503.1
    Z78503.1
    Z78502.1
    Z78502.1
    Z78501.1
    Z78501.1
    Z78500.1
    Z78500.1
    Z78499.1
    Z78499.1
    Z78498.1
    Z78498.1
    Z78497.1
    Z78497.1
    Z78496.1
    Z78496.1
    Z78495.1
    Z78495.1
    Z78494.1
    Z78494.1
    Z78493.1
    Z78493.1
    Z78492.1
    Z78492.1
    Z78491.1
    Z78491.1
    Z78490.1
    Z78490.1
    Z78489.1
    Z78489.1
    Z78488.1
    Z78488.1
    Z78487.1
    Z78487.1
    Z78486.1
    Z78486.1
    Z78485.1
    Z78485.1
    Z78484.1
    Z78484.1
    Z78483.1
    Z78483.1
    Z78482.1
    Z78482.1
    Z78481.1
    Z78481.1
    Z78480.1
    Z78480.1
    Z78479.1
    Z78479.1
    Z78478.1
    Z78478.1
    Z78477.1
    Z78477.1
    Z78476.1
    Z78476.1
    Z78475.1
    Z78475.1
    Z78474.1
    Z78474.1
    Z78473.1
    Z78473.1
    Z78472.1
    Z78472.1
    Z78471.1
    Z78471.1
    Z78470.1
    Z78470.1
    Z78469.1
    Z78469.1
    Z78468.1
    Z78468.1
    Z78467.1
    Z78467.1
    Z78466.1
    Z78466.1
    Z78465.1
    Z78465.1
    Z78464.1
    Z78464.1
    Z78463.1
    Z78463.1
    Z78462.1
    Z78462.1
    Z78461.1
    Z78461.1
    Z78460.1
    Z78460.1
    Z78459.1
    Z78459.1
    Z78458.1
    Z78458.1
    Z78457.1
    Z78457.1
    Z78456.1
    Z78456.1
    Z78455.1
    Z78455.1
    Z78454.1
    Z78454.1
    Z78453.1
    Z78453.1
    Z78452.1
    Z78452.1
    Z78451.1
    Z78451.1
    Z78450.1
    Z78450.1
    Z78449.1
    Z78449.1
    Z78448.1
    Z78448.1
    Z78447.1
    Z78447.1
    Z78446.1
    Z78446.1
    Z78445.1
    Z78445.1
    Z78444.1
    Z78444.1
    Z78443.1
    Z78443.1
    Z78442.1
    Z78442.1
    Z78441.1
    Z78441.1
    Z78440.1
    Z78440.1
    Z78439.1
    Z78439.1



```python
#before these complements can be saved to a file, we first need to make some SeqRecord objects
#this can be done with a generator expression
records = (
    rec.reverse_complement(id="rc_" + rec.id, description="reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta", "fasta")
    if len(rec) < 700 #add a conditional
)
SeqIO.write(records, "rev_comp.fasta", "fasta")
```




    18


## Sequence Alignments
An overview of the Align library.

```python
#import necessary libraries/functions
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
```


```python
#let's start with a single sequence alignment using a stockholm file for Phage_Coat_Gp8 (PF05371)
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/PF05371_seed.sth
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
#as you can see, the sequences are truncated
#we can format it our own way as follows
print("Alignment length %i" % alignment.get_alignment_length())
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    Alignment length 52
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
#the file also includes some database cross-references to the Protein Data Bank and associated secondary structures
#we can access them with an iterator
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
```

    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']



```python
#we can also look at sequence annotations
for record in alignment:
    print(record)
```

    ID: COATB_BPIKE/30-81
    Name: COATB_BPIKE
    Description: COATB_BPIKE/30-81
    Database cross-references: PDB; 1ifl ; 1-52;
    Number of features: 0
    /accession=P03620.1
    /start=30
    /end=81
    Per letter annotation for: secondary_structure
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
    ID: Q9T0Q8_BPIKE/1-52
    Name: Q9T0Q8_BPIKE
    Description: Q9T0Q8_BPIKE/1-52
    Number of features: 0
    /accession=Q9T0Q8.1
    /start=1
    /end=52
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
    ID: COATB_BPI22/32-83
    Name: COATB_BPI22
    Description: COATB_BPI22/32-83
    Number of features: 0
    /accession=P15416.1
    /start=32
    /end=83
    Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
    ID: COATB_BPM13/24-72
    Name: COATB_BPM13
    Description: COATB_BPM13/24-72
    Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
    Number of features: 0
    /accession=P69541.1
    /start=24
    /end=72
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPZJ2/1-49
    Name: COATB_BPZJ2
    Description: COATB_BPZJ2/1-49
    Number of features: 0
    /accession=P03618.1
    /start=1
    /end=49
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
    ID: Q9T0Q9_BPFD/1-49
    Name: Q9T0Q9_BPFD
    Description: Q9T0Q9_BPFD/1-49
    Database cross-references: PDB; 1nh4 A; 1-49;
    Number of features: 0
    /accession=Q9T0Q9.1
    /start=1
    /end=49
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPIF1/22-73
    Name: COATB_BPIF1
    Description: COATB_BPIF1/22-73
    Database cross-references: PDB; 1ifk ; 1-50;
    Number of features: 0
    /accession=P03619.2
    /start=22
    /end=73
    Per letter annotation for: secondary_structure
    Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')



```python
#now that we've read an alignment file, let's try writing one
#we'll start by creating a few MultipleSeqAlignment objects using SeqRecord objects
align1 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"),
        SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"),
        SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma"),
    ]
)

align2 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("GTCAGC-AG"), id="Delta"),
        SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
        SeqRecord(Seq("GTCAGCTAG"), id="Zeta"),
    ]
)

align3 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"),
        SeqRecord(Seq("ACTAGTACAGCT-"), id="Theta"),
        SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"),
    ]
)
my_alignments = [align1, align2, align3]
```


```python
#now we can write them into a PHYLIP file
AlignIO.write(my_alignments, "my_example.phy", "phylip")
```




    3




```python
#let's take a look at it
alignments = AlignIO.parse("my_example.phy", "phylip")
for alignment in alignments:
    print(alignment, "\n")
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma 
    
    Alignment with 3 rows and 9 columns
    GTCAGC-AG Delta
    GACAGCTAG Epsilon
    GTCAGCTAG Zeta 
    
    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota 
    



```python
#we can also look at individual alignments within this file
#first, we need to convert it into a list
alignments = list(AlignIO.parse("my_example.phy", "phylip"))
```


```python
#then we can define and print individual alignments
first_align = alignments[0]
print(first_align)
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma



```python
last_align = alignments[-1]
print(last_align)
```

    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota



```python
#we can convert between sequence alignment file formats in the same way we convert between sequence file formats
#via the convert function
count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
#or with the parse and write functions
alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
#the above examples converted a stockholm file to a clustal file
#let's try converting to PHYLIP instead
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip")
```




    1




```python
#however, it should be noted that the normal PHYLIP format truncates the sequence identifiers at 10 characters which can inhibit readability
#we can get around this by using a more relaxed variant instead
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip-relaxed")
```




    1




```python
#we can manipulate the identifiers, which is useful for if you have to work within the 10 character confines of a normal PHYLIP format sequence
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
name_mapping = {}
for i, record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" % i
print(name_mapping)
```

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}



```python
#the above example appends a "seq" identifier to the beginning of each sequence
AlignIO.write([alignment], "PF05371_seed.phy", "phylip")
```




    1




```python
#we can slice alignment sequence objects in a similar manner to slicing sequence or sequence record objects
#first let's set up our alignment for this example
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
print("Number of rows: %i" % len(alignment))
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    Number of rows: 7
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
#then we can slice it
print(alignment[3:7])
```

    Alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
#we can select my column using a double index, in like manner to a number matrix in NumPy
print(alignment[2, 6])
```

    T



```python
#a shorthand version of the above line could be
print(alignment[2].seq[6])
```

    T



```python
#you can pull a string out of single column
print(alignment[:, 6])
```

    TTT---T



```python
#you can select a range of columns
print(alignment[3:6, :6])
```

    Alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49



```python
#if we leave the first index blank, it will print all the rows
print(alignment[:, :6])
```

    Alignment with 7 rows and 6 columns
    AEPNAA COATB_BPIKE/30-81
    AEPNAA Q9T0Q8_BPIKE/1-52
    DGTSTA COATB_BPI22/32-83
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49
    FAADDA COATB_BPIF1/22-73



```python
#we can use this method to remove a section
print(alignment[:, 6:9])
```

    Alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73



```python
#or slice everything after a given column - in this case, the ninth
print(alignment[:, 9:])
```

    Alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
#you can add alignment objects together as well - it works column by column
edited = alignment[:, :6] + alignment[:, 9:]
print(edited)
```

    Alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
#you can also use this to combine alignments from several different genes, so long as the identifiers line up
#sorting the alignment rows alphabetically helps with this
edited.sort()
print(edited)
```

    Alignment with 7 rows and 49 columns
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49



```python
#the substitutions property reports how often letters in the alignment are substituted for one another
alignment = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTCCTA"), id="seq1"),
        SeqRecord(Seq("AAT-CTA"), id="seq2"),
        SeqRecord(Seq("CCTACT-"), id="seq3"),
        SeqRecord(Seq("TCTCCTC"), id="seq4"),
    ]
)
print(alignment)
```

    Alignment with 4 rows and 7 columns
    ACTCCTA seq1
    AAT-CTA seq2
    CCTACT- seq3
    TCTCCTC seq4



```python
substitutions = alignment.substitutions
print(substitutions)
```

        A    C    T
    A 2.0  4.5  1.0
    C 4.5 10.0  0.5
    T 1.0  0.5 12.0
    



```python
#you can use the select method to add missing letters and change the order
m = substitutions.select("ATCG")
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
#the actual order of the pairs is completely arbitrary
m = substitutions.select("ACTG")
print(m)
```

        A    C    T   G
    A 2.0  4.5  1.0 0.0
    C 4.5 10.0  0.5 0.0
    T 1.0  0.5 12.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
#ClustalW is a tool for running multiple sequence alignments from the command line
from Bio.Align.Applications import ClustalwCommandline
cline = ClustalwCommandline("clustalw2", infile="opuntia.fasta")
print(cline)
```

    clustalw2 -infile=opuntia.fasta



```python
#let's try it out with these example files
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.aln
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.dnd
align = AlignIO.read("opuntia.aln", "clustal")
print(align)
```

    Alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191



```python
#we can also generate phylogenetic trees
from Bio import Phylo
tree = Phylo.read("opuntia.dnd", "newick")
Phylo.draw_ascii(tree)
```

                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658
    
```python
#import necessary libraries/functions
from Bio import Align
from Bio.Seq import reverse_complement
```


```python
#pairwise sequence alignment aligns sequences by assigning similarity scores to them and using them to find the optimal alignment
aligner = Align.PairwiseAligner()
```


```python
#you can set the alignment parameters in the constructor
aligner = Align.PairwiseAligner(match_score=1.0)
```


```python
#or after the fact
aligner.match_score = 1.0
```


```python
#let's use it to calculate the alignment score between two sequences
target = "GAACT"
query = "GAT"
score = aligner.score(target, query)
score
```




    3.0




```python
#to see the actual alignments themselves
alignments = aligner.align(target, query)
for alignment in alignments:
    print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    
    target            0 GAACT 5
                      0 |-|-| 5
    query             0 G-A-T 3
    



```python
#that was a global alignment
#we can also perform a local alignment
aligner.mode = "local"
target = "AGAACTC"
query = "GAACT"
score = aligner.score(target, query)
score
```




    5.0




```python
alignments = aligner.align(target, query)
for alignment in alignments:
    print(alignment)
```

    target            1 GAACT 6
                      0 ||||| 5
    query             0 GAACT 5
    



```python
#all alignment parameters are stored as part of the PairwiseAlignment object
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: 0.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: local
    



```python
#we can check what algorithm is in use
aligner.algorithm
```




    'Smith-Waterman'




```python
#and the significant difference - any difference less than the epsilon value will be considered insignificant (i.e. the two scores are considered functionally equal to one another)
aligner.epsilon
```




    1e-06




```python
#let's go back to our original parameters
target = "GAACT"
query = "GAT"
alignments = aligner.align(target, query)
alignment = alignments[0]
alignment
```




    <Alignment object (2 rows x 5 columns) at 0x7f12f84ac690>




```python
#we can retrieve values individually
#for example, alignment score
alignment.score
```




    3.0




```python
#alignment target
alignment.target
```




    'GAACT'




```python
#alignment query
alignment.query
```




    'GAT'




```python
#the coordinates
alignment.coordinates
```




    array([[0, 2, 4, 5],
           [0, 2, 2, 3]])




```python
#the length
len(alignment)
```




    2




```python
#or the shape
alignment.shape
```




    (2, 5)




```python
#as previously shown, we can also set the mode
#notice how it clips off the extraneous T's at the end of the matched sequences
aligner.mode = "local"
local_alignments = aligner.align("TGAACT", "GAC")
local_alignment = local_alignments[0]
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
#we can set the mismatch score
aligner = Align.PairwiseAligner(match = 1.0, mismatch_score = -10)
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: -10.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: global
    



```python
alignments = aligner.align("AAACAAA", "AAAGAAA")
len(alignments)
```




    2




```python
#notice how the aligner introduces gaps into the sequences - this is part of its optimization process when finding mismatches according to our set mismatch score
print(alignments[0])
```

    target            0 AAAC-AAA 7
                      0 |||--||| 8
    query             0 AAA-GAAA 7
    



```python
#let's print it
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
#we can sort alignments easily with the sort method
local_alignment.sort()
print(local_alignment)
```

    target            0 GA-C 3
                      0 ||-| 4
    query             1 GAAC 5
    



```python
#we can combine it with the reverse complement method from Seq if needed
target = "AAACCC"
query = "AACC"
aligner = Align.PairwiseAligner(mismatch_score = -1, internal_gap_score = -1)
```


```python
#here's the normal alignment
aligner.score(target, query)
```




    4.0




```python
#here's the reverse complement alignment
aligner.score(target, reverse_complement(query))
```




    0.0




```python
#the scores are reversed if we read it as part of the negative strand
#here's the reverse complement alignment
aligner.score(target, reverse_complement(query), strand = "-")
```




    4.0




```python
#here's the normal alignment
aligner.score(target, query, strand = "-")
```




    0.0




```python
alignments = aligner.align(target, query)
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
#we can export this in the BED format
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	+	1	5	0	1	4,	0,
    



```python
#let's read the regular sequence (not the reverse complement) as a negative strand
alignments = aligner.align(target, query, strand = "-")
len(alignments)
```




    2




```python
#there are no alignments at all, so it can't be exported to BED
print(alignments[0])
```

    target            0 AAACCC----  6
                      0 ---------- 10
    query             4 ------GGTT  0
    



```python
#we can fix this by setting the gap scores
aligner.left_gap_score = -0.5
aligner.right_gap_score = -0.2
aligner.score(target, query)
```




    3.3




```python
alignments = aligner.align(target, query)
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
#let's try it with the reverse complement read as the negative strand
alignments = aligner.align(target, reverse_complement(query), strand = "-")
print(alignments)
```

    <Bio.Align.PairwiseAlignments object at 0x7f12f84aef90>



```python
#the score matches the regular, positive strand
aligner.score(target, reverse_complement(query), strand = "-")
```




    3.3




```python
#the alignments remain the same
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             4 -AACC- 0
    



```python
#predictably, reading the reverse complement as a positive strand gives a worse result
alignments = aligner.align(target, reverse_complement(query), strand = "+")
aligner.score(target, reverse_complement(query), strand = "+")
```




    -2.0


## OpenCV
A look at the Python computer visualization suite, OpenCV.

```python
#import necessary libraries/functions
import numpy as np
import cv2
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
#load the image being used
img = cv2.imread("mushrum.jpg")
```


```python
#we can check its type - it's a NumPy array
type(img)
```




    numpy.ndarray




```python
#note that an incorrect path does not return an error
img_wrong = cv2.imread("wrong/path/doesnot/abcdegh.jpg")
```


```python
#but if we check its type...
type(img_wrong)
```




    NoneType




```python
#we can display the image now
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fe114e4cb50>




![output_5_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/86175bbe-2014-417c-9101-043dd914f4ab)




```python
#the color channels are out of order - most images use an RGB triplet palette in that order, but OpenCV interprets channel values in the order BGR instead
#we'll need to reorder the channels
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7fe11457d5d0>




![output_6_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/01ec1185-7b24-4bc1-9905-c096c570f30b)




```python
#we can display the image's resolution and number of color channels with the shape method
fix_img.shape
```




    (539, 512, 3)




```python
#we can also convert it to a simple monochrome channel
img_gray = cv2.imread("mushrum.jpg", cv2.IMREAD_GRAYSCALE)
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7fe1144e7d50>




![output_8_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/39f8ab18-8907-4dcb-b65c-3717e53c4889)




```python
#again, the first channel - in this case, the only channel - is interpreted as blue
#we can fix this by specifying a colormap
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fe11423f450>




![output_9_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/d7db3c4d-bc08-490a-b11d-592eea8e9cb4)




```python
#the image can be resized
new_img = cv2.resize(fix_img, (1000, 400))
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7fe114220710>



![output_10_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/f710e129-e413-4094-b951-253f2e15ce82)



```python
#we can also specify an aspect ratio instead
w_ratio = 0.5
h_ratio = 0.5
new_img = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7fe1141798d0>




![output_11_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/600e389b-8f52-48fa-b4fd-f683fa880759)




```python
#here we can see that the image has been scaled down by half while still maintaining its original aspect ratio
new_img.shape
```




    (270, 256, 3)




```python
#we can flip the image along the horizontal axis
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7fe114160810>




![output_13_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/96423cd8-8b2d-47a0-8fb8-268df811b81a)




```python
#the vertical axis
flip_img2 = cv2.flip(fix_img, 1)
plt.imshow(flip_img2)
```




    <matplotlib.image.AxesImage at 0x7fe1140c7990>




![output_14_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/a85c2832-42b9-4a22-b251-712a3927306a)




```python
#or along both axes
flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img2)
```




    <matplotlib.image.AxesImage at 0x7fe11402eb10>




![output_15_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/d89020fe-3e53-44ee-80da-8f54bceb89ef)




```python
#and we can export as a new image
#note that it does not save the color correction, only the new orientation
cv2.imwrite('mushrum_fixed_image.jpg', flip_img)
```




    True




```python
#let's mess with colors some more
#first we'll import the original image
img = cv2.imread("mushrum.jpg")
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fe10cfeec10>




![output_17_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/deb5cf3b-af4a-4b8c-be45-5564a5eba9ad)




```python
#and fix the colors
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7fe10cf53d90>




![output_18_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/757cd2e3-8e00-44e2-b79d-087918492525)




```python
#color values are stored in RGB format for this image, but we can convert it to HSV instead
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7fe10cebaed0>




![output_19_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/37472fae-ebb3-488c-b7cc-4bf88e77b386)




```python
#or HSL
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
plt.imshow(img3)
```




    <matplotlib.image.AxesImage at 0x7fe10ceae0d0>




![output_20_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/dbca6a9c-f40b-4cf0-bc00-0400dd047dd7)




```python
#let's look at alpha compositing
#we'll use two images this time
img1 = cv2.imread('deeprock.jpg')
img2 = cv2.imread('mushrum.jpg')
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7fe10c756fd0>




![output_22_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/68830841-2813-43d3-85f8-84824b369d99)



```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7fe10c9a2550>




![output_23_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/b4791d02-8a2c-4ca6-afa8-fb5b1093f52d)




```python
#convert colors first
#this isn't really necessary for the img1 I chose since it's just black and white, but I might as well to follow along with the course
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7fe10ca4fed0>




![output_25_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/0e402df7-7c14-4cff-8257-618216bd1d7e)




```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7fe10c9e2910>




![output_26_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/d86034e7-0b52-4f06-a3cc-e696722c763a)




```python
#we'll set them both to the same resolution
img1 = cv2.resize(img1, (1200, 1200))
img2 = cv2.resize(img2, (1200, 1200))
```


```python
#to set individual transparencies, we adjust the alpha channels for the two component images
#don't let the names confuse you - the alpha component and the beta component both have alpha channels (which encode transparency information)
alpha = 0.5
beta = 0.5
```


```python
#then we can composite the two images together
#no gamma correction is being applied
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma = 0)
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7fe10c2b2e50>




![output_29_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/2d384d2c-951e-49a4-84a9-2a88dbf9ec2c)




```python
#we can play with the transparencies
#alpha channel values are between 0 and 1
alpha = 0.8
beta = 0.2
blended1 = cv2.addWeighted(img1, alpha, img2, beta, gamma = 0)
plt.imshow(blended1)
```




    <matplotlib.image.AxesImage at 0x7fe10c21c250>




![output_30_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/9644634c-fd6c-4a7d-bab5-dcc08b33a386)




```python
#we can adjust the size and position of the individual component images
#first let's reload and correct the original images
img1 = cv2.imread('deeprock.jpg')
img2 = cv2.imread('mushrum.jpg')

img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
#then we'll resize just img1
img1 = cv2.resize(img1, (200,200))
```


```python
#let's relabel the images for readability purposes
large_img = img2
small_img = img1
```


```python
#now we can apply an offset
x_offset = 0
y_offset = 0
x_end = x_offset + small_img.shape[1]
y_end = y_offset + small_img.shape[0]
#and composite
large_img[y_offset:y_end, x_offset:x_end] = small_img
plt.imshow(large_img)
```




    <matplotlib.image.AxesImage at 0x7fe10c20fc90>



![output_34_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/61bc9c9d-8b57-4f10-b580-c41fa6486553)




```python
#now let's look at setting thresholds and processing images
#we'll use crossword.jpg and rainbow.jpg - link below
#https://github.com/worklifesg/Python-for-Computer-Vision-with-OpenCV-and-Deep-Learning/tree/main/3.%20Image%20Processing
```


```python
img = cv2.imread('rainbow.jpg')
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fe10c6956d0>




![output_36_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/f3f8fb96-c203-4a49-a474-5646912438b5)




```python
#let's start with thresholding - for example, cancelling out unwanted background colors
#first we'll desaturate it
img = cv2.imread('rainbow.jpg', 0)
plt.imshow(img, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fe10c677350>



![output_37_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/7d111cdb-188b-44c4-99b2-234a96ad4fb0)




```python
#we set the threshold to 127 since there are a total of 255 color values in the image
#any value below 127 will be one color, whereas any value above 127 will be another
ret1, thresh1 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fe10c5d90d0>




![output_38_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/108a2002-3980-494b-a1d2-41c002654963)




```python
#we can invert these limits - the resulting image is an inverse of the grayscale rainbow above
img2 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fe10c5b3dd0>




![output_39_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/c5380ca0-cdc1-4253-bf7e-a1ba2d1171d0)




```python
#we can combine the first two images to allow for some leniency in the thresholding
img3 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img3, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fe10c485590>




![output_40_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/6996d91c-62d3-45b1-98b8-0f12ab2a813b)




```python
#let's apply this to the crossword puzzle
img_r = cv2.imread('crossword.jpg', 0)
plt.imshow(img_r, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fe10c461e90>




![output_41_1](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/c5bd0670-fd4e-41ae-967a-fc7a8c292a35)




```python
#to make things easier, let's create a function to quickly desaturate and resize the image
def show_pic(img):
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
show_pic(img_r)
```



![output_43_0](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/2e0bd72e-ca62-410e-b3c2-adf7546ed1db)



```python
#we'll apply a binary threshold
ret, th1 = cv2.threshold(img_r, 127, 255, cv2.THRESH_BINARY)
show_pic(th1)
```



![output_44_0](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/e9aaf2db-7877-407d-b62a-165769f29c89)



```python
#adjust the threshold a bit
ret, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```



![output_45_0](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/45581884-2af5-43f3-a959-94c3ea4383b8)



```python
#alternatively
th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
show_pic(th2)
```



![output_46_0](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/9ca6fb52-885d-4803-a411-e56df95ff4c7)



```python
#now we'll layer the two images to get a clearer picture
blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                         src2 = th2, beta = 0.4, gamma = 0)
show_pic(blended)
```



![output_47_0](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/a0adbbbb-cb06-4e27-9a9c-fdaca5780b26)



```python
#alternatively
th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)

blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                         src2 = th3, beta = 0.4, gamma = 0)
show_pic(blended)
```


![output_48_0](https://github.com/rjpellegr/RJP023_Advanced_Python_Portfolio/assets/134185456/1480148d-bf07-4092-969a-0fe93031810d)



```python

```
