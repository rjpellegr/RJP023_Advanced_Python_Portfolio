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
