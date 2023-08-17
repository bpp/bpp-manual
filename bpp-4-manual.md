# BPP Documentation

## BPP Explained

### Overview

Bayesian Phylogenetics & Phylogeography (BPP) is a Bayesian Markov chain Monte Carlo (MCMC) program
for analyzing DNA sequence alignments of multiple loci from multiple
closely-related species under the multispecies coalescent (MSC) model
([Yang 2002](https://doi.org/10.1093/genetics/162.4.1811); [Rannala and Yang 2003](https://doi.org/10.1093/genetics/164.4.1645)). See [Xu 2016](https://doi.org/10.1534/genetics.116.190173) and [Rannala et al 2020](https://inria.hal.science/PGE/hal-02535622) for reviews
of the MSC model. The program requires 3 input files:

1.  Sequence data file (see [Sequence File](#sequence-file)).

2.  Imap file specifying the population source of each sequence
    (see [Imap file](#imap-file)).

3.  Control file specifying other variables needed to run the program
    (e.g., initial species tree, prior distributions,
    etc)(see [Control File](#control-file)).

This manual will provide details for creating these input files
(describing syntax and permissible values for variables). The input
files should be in text format and edited using a text format editor
(for example, a Unix editor such as emacs, or vi). A basic set of map
and control files using sensible default priors based on the sequence
data can also be created using the web tool *Minimalist BPP* available
at

<https://brannala.github.io/bpps/> 

These auto-produced files can
then be used as a template for more editing.

### Possible types of analyses

The program can be used to conduct four different types of analyses,
each specified using two variables in the control file:

-   **A00** (`speciesdelimitation = 0`, `speciestree = 0`): estimation
    of species divergence times and population sizes under the MSC, IM
    and MSC-I models when the species tree is specified by the user
	([Rannala and Yang 2003](https://doi.org/10.1093/genetics/164.4.1645); [Flouri et al 2020](https://inria.hal.science/PGE/hal-02536475))


-   **A01** (`speciesdelimitation = 0`, `speciestree = 1`): inference of
    the species tree when sequences are assigned to species by the user
    ([Rannala and Yang 2017](https://doi.org/10.1093/sysbio/syw119))

-   **A10** (`speciesdelimitation = 1`, `speciestree = 0`): species
    delimitation using a user-specified guide tree
    ([Yang and Rannala 2010](https://doi.org/10.1073/pnas.0913022107); [Rannala and Yang 2013](https://doi.org/10.1534/genetics.112.149039))

-   **A11** (`speciesdelimitation = 1`, `speciestree = 1`): joint
    species delimitation and species tree inference using unguided
    species delimitation [Yang and Rannala 2014](https://doi.org/10.1093/molbev/msu279)

A detailed guide to the use of each of these 4 methods of analysis, as
well as example control files, are provided in
[Methods of Analysis](#methods-of-analysis). The BPP tutorial also
provides examples illustrating the four types of analyses
([Yang 2015](https://doi.org/10.1093/czoolo/61.5.854); [Flouri et al 2020](https://inria.hal.science/PGE/hal-02536475)).

### Model Parameters

The basic parameters in the MSC model include the **species divergence times** 
( $\tau$ ) and mutation scaled **population sizes**

$\theta = 4N\mu$

where $N$ is the effective population size and $\mu$
is the mutation rate per site per generation. Note that $\theta$
specifies the average proportion of sites that have different bases when
comparing two sequences sampled at random from the population. Both
$\tau$ and $\theta$ are measured in units of expected number of
mutations per site ( $\tau$ can be converted to units of years if the
substitution rate per site per year is known). For example, if $\tau$ is
an estimated divergence time in units of expected substitutions per site
(from a BPP analysis) and the substitution rate per site per year is
$\mu$, then the estimated divergence time in years is, $\tau'$, is

$\tau' = \frac{\tau}{\mu}$ 

and the diploid effective population size is

$N = \frac{\theta}{4\mu}$ 

For a species tree with $s$ species,
there are ($s – 1$) divergence times ($\tau$) and at most ($2s – 1$)
population size parameters ($\theta$). Analysis A00 estimates those
parameters when the species delimitation and species tree are fixed
(specified by the user). Analyses A01, A10, and A11 estimate parameters
and model probabilities of the different MSC models. The
multispecies-coalescent-with-introgression (MSC-I) model, first
implemented in BPP version 4.1 and described in
[The MSC-I Model](#the-msc-i-model), adds introgression nodes to the tree, each
with a new parameter, the **introgression probability** ($\varphi$) see [Flouri et al 2020](https://doi.org/10.1093/molbev/msz296). 
The multispecies-coalescent-with--migration (MSC-M) model, first implemented in BPP version 4.6 and described in [The MSC-M Model](#the-msc-m-model), adds a matrix of **instantaneous migration rates** (**M**). 
For reviews of the MSC model, see [Yang 2014](https://global.oup.com/academic/product/molecular-evolution-9780199602605?lang=es&cc=gb) (Chapter 9), [Xu 2016](https://doi.org/10.1534/genetics.116.190173) and [Rannala et al 2020](https://inria.hal.science/PGE/hal-02535622).

### Model Assumptions

The MSC model implemented in BPP makes two basic assumptions about the
data: *no recombination between sites of the same locus*, and *free
recombination between loci*. The original MSC model ([Rannala and Yang 2003](https://doi.org/10.1093/genetics/164.4.1645)) also
assumes complete genetic isolation between populations, but more recent
implementations allow gene flow between species (either continuous
migration or episodic introgression); the MSC-I model implemented in BPP
allows episodic introgression [Flouri et al 2020](https://doi.org/10.1093/molbev/msz296) and the IM model allows
continuous migration (see
[Introgression and Migration Models](#introgression-and-migration-models) and [Isolation With Migration Models](#isolation-with-migration-models)). The default model assumes a strict
molecular clock (constant substitution rates among lineages) and
relaxed-clocks (with substitution rates varying among lineages) are
implemented through variable clock models (see
[Substitution Models](#substitution-models)).

### Nuclear genomic data

Ideal data for analysis using the BPP program are loosely-linked short
genomic segments. Such data are likely to satisfy the model assumptions
because intra-locus recombination is rare for short segments (500 to
1000bp) and distantly spaced segments undergo frequent recombination and
thus have nearly independent histories. Genealogical trees are assumed
to be a product of neutral evolution (not influenced by natural
selection). However, protein-coding gene sequences appear to be useable
in BPP analyses since most proteins are performing similar
functions in closely related species and the main effect of purifying
selection on nonsynonymous mutations is a reduction of the neutral
mutation rate ([Shi and Yang 2018](https://doi.org/10.1093/molbev/msx277); [Thawornwattana et al 2018](https://doi.org/10.1093/molbev/msy158)). 
It is a good idea to separate the noncoding and coding regions of the genome into two
datasets.

The program allows the option of either a constant rate (strict
molecular clock) or variable rates of substitution (relaxed molecular
clocks) among lineages. Different substitution models are available,
with the most complex (parameter rich) being the General Time-Reversible
(GTR) model and the simplest the Jukes-Cantor (JC69) mutation model. All
the models correct for multiple hits at individual sites (see [Substitution Models](#substitution-models)). 
The use of the JC69 model and the strict molecular clock model should be limited to
closely related species with sequence divergence not much higher than
10%. Sequences from distantly related species should be analyzed using a
GTR substitution model in combination with one of the relaxed clock
models (see [Among-species rate variation](#among-species-rate-variation). Models allowing substitution rate variation among sites 
(see [Among-site rate variation](#among-site-rate-variation)) and
among loci (see [Among-locus rate variation](#among-locus-rate-variation))are also implemented in BPP (see [Models of Substitution Rate Variation](#models-of-substitution-rate-variation)). 

### Mitochondrial data

The MSC model accommodates genealogical heterogeneity across the genome.
Thus, different regions of the autosomal genome may have different gene
tree topologies and branch lengths (coalescent times). While the
mutational process may also vary along the genome, this heterogeneity is
believed to be much less important. Thus multiple genes from the
mitochondrial genome, which in most species does not undergo
recombination, should be treated as one 'locus' in the MSC-based
analysis. The mitochondrial genome often has a mutation rate that
differs from the nuclear genome, as well as a different effective
population size from that of autosomes. While the model of locus-rate
variation (option variable `locusrate`) and the heredity scalar (option
variable `heredity`) are designed to deal with this, it may be prudent
to analyze the loci from the nuclear genome separately from the single
mitochondrial locus.

### Identical sequences and phylogenetic signal

The model implemented in BPP assumes that the sequences
represent random samples from the different species. *Sequences from the
same species that are identical should all be used*. It is incorrect to
use only the unique haplotypes, which will lead to biased parameter
estimates. Similarly, it is incorrect to filter loci based on bootstrap
support values and use only those loci with a high "phylogenetic
information" signal.

### How to Cite BPP

If you publish results obtained using BPP you should cite the version
number of the BPP release you used as well as the canonical citation,
which as of release v4.3.8 is
```
Tomáš Flouri, Xiyun Jiao, Bruce Rannala, Ziheng Yang (2018) Species
Tree Inference with BPP Using Genomic Sequences and the Multispecies
Coalescent, Molecular Biology and Evolution 35: 2585--2593.
```
You can also cite the one of the BPP tutorials and the
original papers describing the methods you used (see above). Be sure to
state the priors that you used since they are necessary for
reproducibility. If you conduct a joint analysis of species delimitation
and species tree inference, your method description may look like the
following (replace the specific values in bold with those you used):

> "Joint Bayesian species delimitation and species tree estimation was conducted using the program
> BPP (Flouri et al., 2018; Yang, 2015). The method uses the multispecies coalescent model to compare
> different models of species delimitation (Yang and Rannala, 2010; Rannala and Yang, 2013) and species
> phylogeny (Yang and Rannala, 2014; Rannala and Yang, 2017) in a Bayesian framework, accounting
> for incomplete lineage sorting due to ancestral polymorphism and gene tree-species tree discordance.
> The population size parameters ($\theta$s) are assigned the **inverse-gamma prior IG(3, 0.002)**, with mean
> **0.002/(3–1) = 0.001**. The divergence time at the root of the species tree ($\tau_0$) is assigned the 
> **inverse-gamma prior IG(3, 0.004)**, with mean **0.002**, while the other divergence time parameters are specified
> by the uniform Dirichlet distribution (Yang and Rannala, 2010, eq. 2). Each analysis is run at least twice
> to confirm consistency between runs.”

### Features New To BPP 4.0
- Multiple threads
- Introgression models (MSC-I)
- Migration models (MSC-M)
- Local-clock and locus-rates models (clock and locusrate)
- Topological constraints (constraint and outgroup)

## Installing and Running BPP

This manual applies to BPP versions 4.1 and later. The
program is written in C and executables are available for Windows, Mac
OSX and Linux. Alternatively, the BPP program can be
compiled for Unix (Linux, BSD, etc), Mac OSX, or Windows. If you are
using a Windows, Mac OSX or Linux operating system (and you are not a
programmer who prefers to compile from source) you should download the
latest release executables for your machine (go to [Obtaining BPP](#obtaining-bpp)). If you
have not used a command line program before, read
[Using the Command Line](#using-the-command-line) then go to
[Obtaining BPP](#obtaining-bpp). If
you want to compile the program yourself (you must do this if you are
using a Unix variant other than Linux or Mac OSX) go to
[Compiling the BPP Program](#compiling-the-bpp-program). 

### Obtaining BPP

Executables and source code for the most recent BPP
release are always available at:

<https://github.com/bpp/bpp/releases/latest>

Download the distribution file containing an executable for your
operating system and uncompress the file. Change the current directory
to be the root of the subdirectory created by uncompressing the file.
For example, on a Linux machine you could type the following to download
and install bpp version 4.6.2 (the latest version may be different from
this):
```
wget https://github.com/bpp/bpp/releases/download/v4.6.2/bpp-4.6.2-linux-x86_64.tar.gz
tar -xzf bpp-4.6.2-linux-x86_64.tar.gz
cd bpp-4.6.2-linux-x86_64
```
You are now ready to proceed to [BPP Trial Run](#bpp-trial-run)

### Using the Command Line

BPP is a command-line program, so the preferred way of
running it is to open a terminal application and execute the program
from the terminal command line, rather than double-clicking on the
executable file in your file explorer. If you have not used the command
line before, please work through one of the following short tutorials
first:

Windows:  
<http://abacus.gene.ucl.ac.uk/software/CommandLine.Windows.pdf>

Mac OSX:  
<http://abacus.gene.ucl.ac.uk/software/CommandLine.MACosx.pdf>

### Compiling the BPP Program

#### Requirements

To compile the program you will need to have a C compiler and the Make
program installed on your machine. To test whether this software is
installed (on a Unix machine) you can type:
```
cc --version; make --version
```
which should produce output similar to the following if a compiler is
installed:
```{ .yaml .no-copy }
cc (Ubuntu 9.3.0-10ubuntu2) 9.3.0
Copyright (C) 2019 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

GNU Make 4.2.1
Built for x86_64-pc-linux-gnu
Copyright (C) 1988-2016 Free Software Foundation, Inc.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
```
If instead you see an error message such as:
```{ .yaml .no-copy }
Command 'cc' not found
```
then you (or your system) administrator will need to install the
compiler software (using a package manager such as apt on Ubuntu, for
example) before you can compile BPP.

#### Downloading and compiling the source code

If you have the necessary software for compiling you can download the C
source code for the latest release from:

<https://github.com/bpp/bpp/releases/latest>

The link to the source file on github is named `Source code`. You will
need to download and uncompress the file, then change to the source code
subdirectory before executing the commands outlined below. The program
needs to be compiled only once. For example, the following commands use
the gcc compiler to compile the program and move the generated
executable file (bpp) into the bin/ folder.
```
cd bpp
mkdir bin
cd src
make
mv bpp ../bin
```
If you use git you can instead clone the bpp repository and check out
the master branch (which contains source code for the latest stable
version of bpp) and then compile the program:
```
git clone https://github.com/bpp/bpp.git
cd bpp
git checkout master
cd src
make
```
You are now ready to proceed to [BPP Trial Run](#bpp-trial-run)

### BPP Trial Run

Run the program from the command line (rather than double-clicking the
executable) so that you will see any potential error messages. Change
your current directory to the top level of the directory created by
uncompressing the bpp distribution file. For example, in Linux the
following commands will uncompress the distribution file and move you to
the top level of the bpp directory:
```
tar -xvzf bpp-4.6.2-linux-x86_64.tar.gz
cd bpp-4.6.2-linux-x86_64
```
In the bpp/ subdirectory, run the program in Windows by typing the
following command within the Terminal application:
```
bin\bpp --cfile examples\frogs\A00.bpp.ctl
```
In Linux or Mac OSX type the following commands within a terminal:
```
cd examples/frogs/
../../bin/bpp --cfile examples/frogs/A00.bpp.ctl
```
If the program executed successfully, you should see initial output
similar to the following:
```{ .yaml .no-copy }
bpp v4.6.2_linux_x86_64, 31GB RAM, 12 cores
https://github.com/bpp/bpp
	
Auto-selected SIMD ISA: AVX2
	

Starting timer..
Using seed: -1
Parsing species tree... Done
Parsing phylip file... Done
```
Alternatively, you may see one or more error messages. This may be due
to a spelling error, or you may not be in the correct subdirectory. If
this does not appear to be the case and the problem persists you can ask
for help on the BPP discussion group at:

<https://groups.google.com/forum/#!forum/bpp-discussion-group>

When posting on the forum, please specify the exact command you typed
and the error message you received (preferably by copying and pasting
this information from your terminal or attaching a screen shot).

## Input File Formats

Two input files are required for every BPP analysis, the sequence file
and the control file. If more than one species is being analyzed an Imap
file is also required. A few additional input files are required only
when specific options are specified in the control file. Here
we describe the content and format specifications of the sequence, Imap
and control files that will be needed for most BPP analyses.

### Sequence File

The sequence alignments must be in the phylip/paml format, with one
alignment following the other, all in one file. An alignment is called a
locus. Every locus must have at least 2 sequences, different loci can
have different numbers of sequences, and some species may have no
sequences present at some loci. In the phylip/paml format the first line
of the entry for each locus contains two numbers specifying the number
of sequences and the number of sites, respectively, that are present at
the locus. A sequence name must precede each sequence entry and be
separated from it by either a line break or two or more spaces. Each
sequence name must contain an individual ID, which is unique among the
sequences listed for a given locus, preceded by a caret `^` symbol. For
example, the sequence name `GI01234567^Specimen1` has the individual ID
`Specimen1`. The label `GI01234567` that precedes the `^` is optional, so
`^Specimen1` is also a valid sequence name but `Specimen1` is not. The
sequences may be in interleaved format.

DNA bases can be specified using either upper or lower case letters,
missing data can be specified using a `?` symbol, and IUPAC ambiguity
codes may be used. The interpretation of an ambiguity code (as either an
uncertain base or a diploid genotype) depends on whether the data are
specified as phased in the control file. By default (control file
variable `cleandata = 0`, see [Control File](#control-file)),
alignment gaps and ambiguity nucleotides are used in the likelihood
calculation, with gaps treated as question marks 
(see [Yang 2014](https://global.oup.com/academic/product/molecular-evolution-9780199602605?lang=es&cc=gb) pp. 111-112). 
If `cleandata = 1`, all columns with gaps or ambiguity
characters are removed before analysis. An example of a sequence data
file for 2 loci and 4 sequences per locus is given below (you can also
examine the sequence files `frogs.txt` and `yu2001.txt` in the
`examples` subdirectory).
```
    4 20

    dog^alfy1    attgtgccctctctctctca
    dog^alfy2    attgtgccctctctctctca
    cat^smoot    attgtgccctctagctctca
    human^ben    attgtgccctctgactctcc

    4 10

    dog^alfy1    attgtgccct
    dog^alfy2    attgtgccct
    fox^sample1  attgtgccct
    ^ben         attgtgaagt
```
Imap File
---------

Each sequence name has an individual ID tag after the caret `^` symbol.
For example, the sequence name `GI01234567^Specimen1` has the individual
ID `Specimen1`. An **Imap file** maps the individual IDs (specimens) to
their species/populations. For example, the following Imap file maps
`Specimen1` to species A, `Specimen4` to species B, and so on.
```
    Specimen1   A
    Specimen2   A
    Specimen3   A
    Specimen4   B
    Specimen5   C
    Specimen6   C
```
See also the map files `frogs.Imap.txt` and
`Rokas2003-5species-Imap.txt` in the `examples` subdirectory. The MSC
model implemented in BPP infers the relationships among
species/populations not individuals and so it uses only the population
ID for a sequence (such as A, B, ... in the example), and does not use
the individual ID. Also when BPP reads the sequence names,
it use the individual ID to retrieve the species ID for each sequence
but then ignores the sequence name. One motivation for the use of the
Imap file and this two-layer design is that one may wish to analyze the
same sequence data with different population/species assignments of
individual IDs. This can be achieved by editing the small Imap file
rather than editing the larger sequence data file.

![imap](https://github.com/bpp/bpp-manual/blob/main/figures/fig-imap.png?raw=true)  
An example Imap file showing the sequence ID in the Imap file
corresponding to that in the sequence file. Note that the caret `^`
symbol is not allowed in the sequence ID of the Imap file but is
required in the sequence label in the sequence data
file.

Outgroups and Constraints File
------------------------------

The **outgroups and constraints file** is optional and its name is
specified in the control file by setting variable `constraintfile` to
the name of the constraint file. Topological constraints on species
trees can only be specified for analyses A01 (species tree inference
with a fixed species delimitation) and A11 (joint species tree inference
and species delimitation). See [Methods of Analysis](#methods-of-analysis)
for a detailed description of the A01 and
A11 methods of analysis. Note that species trees are always rooted in
BPP, whether you use the clock or relaxed clock. The
control variable `constraintfile` in the control file has the format
```
    constraintfile = constraints.txt
```
Inside the constraint file, three keywords are allowed: `define`,
`constraint`, and `outgroup`. The `define` keyword is used to assign a
name or alias to a clade, `constraint` defines a clade or subtree, and
`outgroup` means that species not on the list (the ingroup species) form
a clade.
```
    define g1 as (G,H);
    # identical to constraint (D,E,F,(G,H));
    constraint (D,E,F,g1);   
    outgroup A,B,C;
    constraint (A,B);
```
Either the equal sign (=) or spaces can be used to separate the keyword
from the specified value. In the example above, there are eight species
in the dataset: A, B, C, D, E, F, G, and H. Line 1 defines an alias `g1`
for clade (G,H), which may be used in subsequent `define` or
`constraint` statements. Line 3 specifies A, B, C as the outgroups; this
is interpreted to mean that species not on the list, (that is, D, E, F,
G and H) form a clade. Since this information is already in line 2, line
3 has no effect in this case.

The modifier keyword `NOT` (case-insensitive) can be used with `define`
to instead specify the complement of the species on the list. Suppose in
the data we have 10 species of ratites (flightless birds), and two
outgroup species, chicken and ostrich. You can specify the constraint as
follows:
```
    define ratites as NOT (chicken, ostrich);
    constraint = (chicken, (ostrich, ratites));
```
Note that BPP interprets `outgroup` to mean that all
species not on the list form a clade. In theory the `outgroup` keyword
is unnecessary because you can always achieve the same thing by using
`define` (especially with the `NOT` option) and `constraint`. However
`outgroup` may be more convenient or explicit in some situations. If
there are four species (A,B,C,D),
```
    outgroup = B,C,D;
```
achieves nothing since all the 15 rooted species trees are allowed,
whereas
```
    outgroup = C,D;
```
means (A,B) form a clade so that 3 rooted species trees are allowed.
Finally if multiple compatible constraints are specified,
BPP merges them into one so that
```
    constraint (E,F,G,H);
    constraint (F,G,H);
    constraint (G,H);
```
is equivalent to
```
    constraint (E,(F,(G,H)));
```
If any two constraints are in conflict, BPP will abort
with an error message. For example the following will cause an error.
```
    constraint = (E,F,G,H);
    constraint = (G,H);
    constraint = (F,G);
```

Control File
------------

The **control file** specifies most aspects of a BPP analysis, including
the names of input and output files, choices of models, prior
distributions on parameters, etc. The format is line oriented, with each
line specifying the value of a particular variable. The default control
file name is `bpp.ctl`. Lines beginning with `*` or `#` are comments and
are ignored by the program. Most often the order of the lines is
unimportant, but there are exceptions (as explained below). Each
non-commented line of the control file contains a variable name and an
assignment of one or more values to the variable. The general syntax of
a line in the control file is:
```
    variable = value
```
where `variable` is the name of a valid variable and `value` is a list
of one or more permissible values for the variable. Most variable
assignments follow this format --- the exceptions are described later.
BPP has complicated dependencies among variables, both in terms of
whether the variable needs to be defined (or has a default value when
undefined) and what its permissible values are. Here we begin by
defining all the variable names, the syntax of their arguments and
permissible types (integers, floating point numbers, etc), and outline
the structure of dependencies among variables

| Symbol/Rule | Definition                                                   |
|-------------|--------------------------------------------------------------|
| b           | Boolean (0,1)                                                |
| +d          | Positive integer                                             |
| f           | floating point number                                        |
| s           | String (or character)                                        |
| t           | Tree in Newick format                                        |
| x*          | List of one or more elements of type x                       |
| y [z]       | z is a set of variables whose form is determined by y        |
| (w,z)       | w and z are two different permissible types                  |
| y [(w,z)]   | whether w or z is the permissible type depends on value of y |

**Table 1.** Symbols used to represent syntax and value types of variables in variable definitions.


Detailed definitions of all the control file variables are
provided at the end of this section, ordered according to their indexes
in Table 2 below. If you are impatient to get started using
the program with worked examples, you might want to now skip ahead to
[Example control file](#example-control-file) which explains the control file
variables in the context of a simple example dataset analysis. You can
later consult this section regarding particular control file variables
as the need arises.

| Index                         | Variable              | Values              | Condition          | Dependencies |
|-------------------------------|-----------------------|---------------------|--------------------|--------------|
| [1](#1-seed)                  | *seed*                | (-d,+d)             | True               | None         |
| [2](#2-usedata)               | *usedata*             | b                   | True               | None         |
| [3](#3-outfile)               | **outfile**           | s                   | True               | None         |
| [4](#4-mcmcfile)              | **mcmcfile**          | s                   | True               | None         |
| [5](#5-seqfile)               | **seqfile**           | s                   | True               | <-17         |
| [6](#6-finetune)              | **finetune**          | b [f*]              | True               | ->8          |
| [7](#7-print)                 | **print**             | b b b b b           | True               | None         |
| [8](#8-burnin)                | `burnin`              | +d                  | True               | <-6          |
| [9](#9-sampfreq)              | **sampfreq**          | +d                  | True               | None         |
| [10](#10-nsample)             | **nsample**           | +d                  | True               | None         |
| [11](#11-species&tree)        | **species&tree**      | +d s* +d* t         | True               | ->{12,16,25} |
| [12](#12-imapfile)            | `Imapfile`            | s                   | 11[1]>1            | <-11         |
| [13](#13-speciesdelimitation) | *speciesdelimitation* | b [b f*]            | True               | ->15         |
| [14](#14-speciestree)         | *speciestree*         | b                   | True               | ->14         |
| [15](#15-speciesmodelprior)   | `speciesmodelprior`   | +d                  | 13[1]=1 OR 14[1]=1 | <-\{13,14\}  |
| [16](#16-phase)               | **phase**             | b*                  | 11[1]              | <-11         |
| [17](#17-nloci)               | *nloci*               | +d                  | True               | ->5 <-31     |
| [18](#18-model)               | *model*               | s [s]               | True               | None         |
| [19](#19-qrates)              | *Qrates*              | b +f +f +f +f +f +f | 18=''7''           | <-18         |
| [20](#20-basefreqs)           | *basefreqs*           | b +f +f +f +f       | 18=''7''           | <-18         |
| [21](#21-alphaprior)          | *alphaprior*          | f f +d              | True               | None         |
| [22](#22-cleandata)           | *cleandata*           | b                   | True               | None         |
| [23](#thetaprior)             | *thetaprior*          | s [(f*, f f s)]     | True               | None         |
| [24](#tauprior)               | *tauprior*            | s f f               | True               | None         |
| [25](#phiprior)               | `phiprior`            | f f                 | 11[4]              | <-11         |
| [26](#locusrate)              | *locusrate*           | +d [(f f f s, s, )] | True               | <-27 ->27    |
| [27](#clock)                  | *clock*               | +d [f f f s s]      | True               | <-26 ->26    |
| [28](#heredity)               | *heredity*            | +d [(f f, s)]       | True               | None         |
| [29](#checkpoint)             | *checkpoint*          | +d [(+d)]           | True               | None         |
| [30](#constraintfile)         | *constraintfile*      | s                   | True               | None         |
| [31](#threads)                | `threads`             | +d [ +d +d ]        | True               | ->17         |

**Table 2.** Complete list of BPP control file variables. The column `Variable` contains the
  control file variable name, the column `Values` specifies (using symbols) the format of the values (symbols are defined in
  Table 1). The column `Condition` describes the conditions (if any) needed for the variable to be required. An entry
  of True in this column means that the variable is always required/defined although possibly not explicitly (there may be a
  default value). The column `Dependencies` lists the variables that either influence the variable (<-) or are influenced by the variable
  (->) with the variables denoted using their index from column `Index`.


Table 2 lists all the variables that can be defined
in the BPP control file. The first column lists an index used to
identify the variable in other columns of the table and to determine the
order in which the variables are defined in subsequent text. The second
column lists the variable name. The variable names are highlighted to
indicate the structure of their dependencies among variables as follows:

-   Bold = variable always needs to be defined, has no default value, and
    does not depend on other variables.

-   Italic = variable always needs to be defined, has a default value if
    no value is provided in control file, and does not depend on other
    variables.

-   Solid background  = whether variable needs to defined (and possibly also its
    permissible values) depends on the value of another variable.

-   Regular = variable always needs to be defined and has no default value
    but its permissible values depend on the value of another variable.

Column 3 lists the permissible types of each variable and the syntax for
specifying arguments in defining the variable. The variable types are
defined in Table 1.

Column 4 lists the condition (state of another variable) that either
determines whether the variable on that line exists (solid background variables)
or determines its permissible values (italic or regular variables). x [y]
refers to value of argument y of the variable with index x. Column 5
lists indexes of variables that either are influenced by the variable
(the notation ->[c,d] indicates that the variable influences
variables c and d) or variables that influence the variable (the
notation <-[a] indicates that variable a influences the variable). If
you are writing a control file from scratch, these dependencies
necessitate that you specify the bold variables first, followed by the
italic (only needed if changing defaults) and then the solid background and regular.


## BPP control file variables

### 1 seed
----------------
```
seed = (-d,+d)
```
**DESCRIPTION**  
Specifies the seed used for the random number generator.  
**VALUES**  
`-d`, use the wall clock (in Windows) or white noise (in Linux) to
generate a seed (store seed in file `SeedUsed`).  
`+d`, use +d as the seed.  
**DEFAULT**  
`-1`  
**COMMENTS**  
Using the same positive integer seed will produce identical results in
different runs, which is useful for debugging, and using different
positive integers in different runs produces different results. Using
$-d$, for example, $-1$, the computer's entropy generator is used as a
source for the seed, and different runs will produce different results.
When evaluating MCMC convergence by comparing results across different
runs be sure to use -d. The positive integer seed automatically
generated when using option -d is stored in a file named SeedUsed and
can be used explicitly to replicate a result. It is recommended that you
run each analysis at least twice using different seeds to confirm that
the results are stable across runs.  
**EXAMPLES**
```
seed = -1
seed = 278
```

### 2 usedata
------------------------------------------------------------------------
```
usedata = b
```
**DESCRIPTION**  
Specifies whether data (likelihood+priors) are used in calculating
probabilities during MCMC or only priors.  
**VALUES**  
`0`, use only the priors to calculate probabilities (likelihood
constant).  
`1`, use likelihood and prior probabilities.  
**DEFAULT**  
`1`  
**COMMENTS**  
The `0` option can be used for debugging, or examining priors. When
using option `0` a MCMC run produces samples from the prior for each
variable.  
**EXAMPLES**
```
usedata = 0
usedata = 1
```

### 3 outfile
------------------------------------------------------------------------
```
outfile = s
```
**DESCRIPTION**  
Sets the path/name of the output file to be the string `s`  
**VALUES**  
`s`, a string of characters specifying the directory path and name
of file that will receive the program output.  
**EXAMPLES**
```
outfile = /home/mickey/mouse_out.txt
```

### 4 mcmcfile
------------------------------------------------------------------------
```
mcmcfile = s
```
**DESCRIPTION**  
Sets the name of the file receiving MCMC output to be the string
`s`  
**VALUES**  
`s`, a string of characters specifying the directory path and name
of file that will receive the MCMC output.  
**EXAMPLES**
```
mcmcfile = /home/mickey/mouse_mcmc.txt
mcmcfile = mcmc.txt
```

### 5 seqfile
------------------------------------------------------------------------
```
seqfile = s
```
**DESCRIPTION**  
Sets the name of the file containing the sequence data to be the string
`s`  
**VALUES**  
`s`, a string of characters specifying the directory path and name
of file that contains the sequence data  
**EXAMPLES**
```
seqfile = /home/mickey/mouse_seq.txt
seqfile = sequences.txt
```

### 6 finetune
------------------------------------------------------------------------
```
finetune = b: [f*]
```
**DESCRIPTION**  
Determines whether step lengths in MCMC proposals are automatically
optimized or manually set to fixed values as well as setting either the
fixed step lengths (manual) or the starting step lengths (optimized).
There are up to 15 proposals with step lengths that can be adjusted. Not
all proposals are defined for any given choice of model and analysis,
however 15 values should always be specified when using manual as the
order is fixed, step lengths specified for unused proposals will be
ignored by the program. The sizes of proposal steps for each parameter
are provided (from left to right) in the same order as listed (from top
to bottom in table. The proposals are as follows:

| Abbreviation | Parameter Proposal Description                                |
|--------------|---------------------------------------------------------------|
| Gage         | Node age on gene tree                                         |
| Gspr         | Subtree pruning and regrafting move                           |
| thet         | Theta                                                         |
| tau          | Tau                                                           |
| mix          | Mixing step (jointly changing mu, tau and gene tree node ages |
| lrht         | empty                                                         |
| phi          | Introgression probability                                     |
| pi           | Stationary frequencies of substitution model                  |
| qmat         | Q matrix elements of substitution model                       |
| alfa         | Alpha parameter of gamma model for rate variation among sites |
| mubr         | empty                                                         |
| nubr         | empty                                                         |
| mu\_i        | empty                                                         |
| nu\_i        | empty                                                         |
| brte         | empty                                                         |


**VALUES**  
`0: Gage Gspr thet tau mix lrht phi pi qmat alfa mubr nubr mu_i nu_i brte`  
Use fixed step lengths in MCMC proposals as specified, where Gage is the
specified step length for proposal Gage, and so on.  
`1: Gage Gspr thet tau mix lrht phi pi qmat alfa mubr nubr mu_i nu_i brte`  
Automatically optimize step lengths in MCMC proposals using specified
initial step lengths.  
**DEFAULT**  
`0`  
Manually fixes all step lengths to be $0.1$.  
`1`  
Fixes all initial proposal step lengths to $0.1$ prior to optimization.  
**DEPENDENCIES**  
If `b` is set to 1 (automatic optimization) then `burnin` has to
be $>200$.  
**COMMENTS**  
The first value, before the colon, is a switch, with 0 meaning no
automatic adjustments by the program and 1 meaning automatic adjustments
by the program. Following the colon are the step lengths for the
proposals used in the program. and then the step lengths specified here
will be the initial step lengths, and the program will try to adjust
them using the information collected during the burnin step. This option
appears to work fine. Some notes about manually adjusting those finetune
step lengths are provided below in section 3.2.  
**EXAMPLES**
```
finetune = 0: .01 .02 .03 .04 .05 .01 .01 .01 .01 .01 .01 .01 .01 .01 .01
finetune = 1
```

### 7 print
------------------------------------------------------------------------
```
print = b b b b b
```
**DESCRIPTION**  
Specifies what information is saved to the output files during the
MCMC.  
**VALUES**  
There are 5 boolean (0,1) variables that specify whether particular
variables are saved to file (1) or not saved (0). From left to right:
variable 1 specifies MCMC samples; variable 2 specifies locus-specific
rate parameters (rate $\mu_i$, variance parameter $\nu_i$ and
species-tree branch rates for locus $r_{ij}$); variable 3 specifies
locus-specific heredity scalars if those are estimated from the data;
variable 4 specifies locus-specific gene trees; and variable 5 specifies
locus-specific substitution-rate parameters (qmat for the $Q$ matrix,
freqs for base frequencies, and alpha for gamma rates for sites).  
**EXAMPLES**
```
print = 0 1 1 1 0
print = 1 1 1 1 1
```

### 8 burnin
------------------------------------------------------------------------
```
burnin = +d
```
**DESCRIPTION**  
Specifies the number of burn-in iterations that are executed in the MCMC
run before sampling begins.  
**VALUES**  
`+d`, a positive integer specifying the number of burn-in
iterations.  
**DEPENDENCIES**  
If `finetune = 1` (automatic optimization) then `burnin` has to
be $>200$.  
**COMMENTS**  
The total number of MCMC iterations is `burnin` + `nsample` $\times$ `sampfreq`.  
**EXAMPLES**
```
burnin = 100
burnin = 10000
```

### 9 sampfreq
------------------------------------------------------------------------
```
sampfreq = +d
```
**DESCRIPTION**  
Specifies the interval at which samples from the MCMC are to be written
to the output file specified by `mcmcfile`.  
**VALUES**  
`+d`, a positive integer specifying the interval between samples
from the MCMC.  
**COMMENTS**  
The total number of MCMC iterations is `burnin` + `nsample` $\times$ `sampfreq`.  
**EXAMPLES**
```
sampfreq = 2
sampfreq = 10
```

### 10 nsample
------------------------------------------------------------------------
```
nsample = +d
```
**DESCRIPTION**  
Specifies the number of samples from the MCMC that are to be written to
the output file specified by `mcmcfile`.  
**VALUES**  
`+d`, a positive integer specifying the number of MCMC samples.  
**COMMENTS**  
The total number of MCMC iterations is `burnin` + `nsample` $\times$ `sampfreq`.  
**EXAMPLES**
```
nsample = 5000
nsample = 25499
```

### 11 species&tree
------------------------------------------------------------------------
```
species&tree = +d s*
+d*
t
```
**DESCRIPTION**  
Specifies the number of species, the species names, the maximum number
of sequences (at any locus) for each species, and the species tree in
Newick format (or species graph in extended Newick format in models with
introgression).  
**VALUES**  
`+d`, a positive integer specifying the total number of species.
`s*`, a list of `+d` strings specifying the name of each
species. `+d*`, a list of `+d` positive integers specifying the
maximum number of sequences (at any locus) for each species. `t`, a
tree (or graph) in Newick format (or Newick extended format) specifying
the relationships among species (and possibly introgression events).  
**DEPENDENCIES**  
The variable `Imapfile` must be specified only when the number of
species `+d` > 1. The length of the boolean array
specified for variable `phase` is equal to the number of species
`+d`. The variable `phiprior` must be specified only if `t`
is an introgression graph.  
**COMMENTS**  
See the subsection *Notation for trees and introgression graphs* for
details on the format of `t`. The role of `t` in the MCMC
analysis depends on the value of the `speciestree` variable. If
`speciestree=1` the variable `t` specifies the starting tree for
the MCMC analysis, otherwise if `speciestree=0` it is the fixed tree
used for estimating parameters or the fixed guide tree for species
delimitation.  
**EXAMPLES**
```
species&tree = 4 H C G O
1 1 1 1
(((H, C), G), O);
```

### 12 Imapfile
------------------------------------------------------------------------
```
Imapfile = s
```
**DESCRIPTION**  
Sets the path/name of the Imap file to be the string `s`.  
**VALUES**  
`s`, a string of characters specifying the directory path and name
of a file that contains the species/sequence map information.  
**DEPENDENCIES**  
The variable `Imapfile` must be specified only when the number of
species defined by the `species&tree` variable
`+d` > 1.  
**COMMENTS**  
See [Imap File](#imap-file) for a complete description of the Imap file format.  
**EXAMPLES**
```
Imapfile = /home/foo/bar_Imap.txt
Imapfile = Imap.txt
```

### 13 speciesdelimitation
------------------------------------------------------------------------
```
speciesdelimitation = b [b f*]
```
**DESCRIPTION**  
Specifies whether species delimitation is performed and if so specifies
the parameters of the species delimitation rjMCMC proposal algorithm.  
**VALUES**  
`0`, specifies no species delimitation.  
`1 0 epsilon`, specifies species delimitation with rjMCMC algorithm
0, with parameter $\epsilon$ as specified in equations 3 and 4 of
[Yang and Rannala 2010](https://doi.org/10.1073/pnas.0913022107). Reasonable values for $\epsilon$ are 1, 2, 5, etc.  
`1 1 alpha m`, specifies rjMCMC algorithm 1, with parameters
$\alpha$ and $m$ as specified in equations 6 and 7 of [Yang and Rannala 2010](https://doi.org/10.1073/pnas.0913022107).
Reasonable values are $\alpha = 1, 1.5, 2$, etc. and $m = 0.5, 1, 2$.  
**DEFAULT**  
`0`  
**DEPENDENCIES**  
If `speciesdelimitation = 1` then `speciesmodelprior` must be defined.  
**EXAMPLES**
```
speciesdelimitation = 1 0 2
speciesdelimitation = 1 1 2 1
```

### 14 speciestree
------------------------------------------------------------------------
```
speciestree = b
```
**DESCRIPTION**  
Specifies whether the species tree is estimated or fixed in the MCMC
analysis.  
**VALUES**  
`0`, specifies that species tree is fixed.  
`1`, specifies that species tree is estimated.  
**DEFAULT**  
`0`  
**DEPENDENCIES**  
If `speciestree = 1` then `speciesmodelprior` must be defined.  
  
**COMMENTS**  
This option invokes the nearest-neighbor interchange (NNI) or
subtree-pruning-regrafting (SPR) algorithm to propose changes to the
species tree topology during the MCMC run.  
**EXAMPLES**
```
speciestree = 1
speciestree = 0
```

### 15 speciesmodelprior
------------------------------------------------------------------------
```
speciesmodelprior
```
**DESCRIPTION**  
Specifies the prior distribution for rooted species trees and species
delimitations.  
**VALUES**  
`0`, specifies equal probabilities for the labeled histories (rooted
trees with the internal nodes ordered by age).  
`1`, specifies equal probabilities for the rooted species trees.  
`2`, specifies equal probabilities for the numbers of species ($1/s$
each for $1,2,\ldots,s$ species given $s$ populations) and divides the
probability for any specific number of species among the compatible
models (of species delimitation and species phylogeny) in proportion to
the number of labeled histories.  
`3`, specifies equal probabilities for the numbers of species ($1/s$
each for $1,2,\ldots,s$ species given $s$ populations) and divides the
probability for any specific number of species among the compatible
models (of species delimitation and species phylogeny) uniformly.  
**DEPENDENCIES**  
The variable `speciesmodelprior` must be defined if either (or both)
`speciestree = 1` or `speciesdelimitation = 1`. If
`speciestree = 1` and `speciesdelimitation = 0` then only
options `speciesmodelprior = 0` or `speciesmodelprior = 1` may
be used.  
**COMMENTS**  
The priors specified by options 2 and 3 are described in [Yang and Rannala 2014](https://doi.org/10.1093/molbev/msu279).
The prior specified by option 3 may be suitable when there are many populations.  
**EXAMPLES**
```
speciesmodelprior = 0
speciesmodelprior = 3
```

### 16 phase
------------------------------------------------------------------------
```
phase = b*
```
**DESCRIPTION**  
Specifies whether the sequence data for each species is phased.  
**VALUES**  
`1`, in position $i$ of the list indicates that data is unphased for
species $i$ in the list of species in variable `species&tree`,
otherwise `0` indicates it is phased.  
**DEFAULT**  
`0`  
**DEPENDENCIES**  
The number of boolean (0,1) variables in the list must equal the number
of species `d+` specified in the `species&tree` variable.  
**COMMENTS**  
If at position $i$ the `phase` variable is 1, each sequence from
species $i$ is treated as an unphased diploid sequence, with
heterozygous sites represented using the ambiguity characters `YRMKSW`.
The missing data tokens `N-?` are allowed and are all treated as ?, but
other ambiguity characters (`HBVD`) are not allowed.

For phasing, each sequence is expanded/resolved into two sequences, and
the program averages over possible phase resolutions according to the
approach of @Gronau2011. For example, a sequence `R...Y` with two
heterozygous sites, `R` (meaning `A/G`) and `Y` (meaning `C/T`), has two
possible phase resolutions:

1.  `A...T` and `G...C`

2.  `A...C` and `G...T`

A sequence with no heterozygous sites will be phased/resolved into two
identical sequences. The MSC model effectively averages over different
phase resolutions of the heterozygous sites in performing the likelihood
calculation. The option may increase memory usage and CPU time
considerably if many sequences exist with many heterozygote sites at a
locus. If the `phase` variable for a species is 0, all sequences
from that species are treated as resolved haplotype sequences, and
ambiguities are interpreted in the usual way. For example, an ambiguity
code `Y` at a position is interpreted to mean that one sequence exists
with an uncertain nucleotide at that position that is either a T or a C.
The program does not allow some sequences from a species to be phased
and other sequences from that same species to be unphased.  
**EXAMPLES**
```
phase = 0 0 0 0 1
phase = 1 1 1 1 1
```

### 17 nloci
------------------------------------------------------------------------
```
nloci = +d
```
**DESCRIPTION**  
Specifies the number of loci present in the sequence data file specified
by variable `seqfile` that will be analysed.  
**VALUES**  
`+d`, a positive integer specifying the number of loci to be
analyzed.  
**DEFAULT**  
If `nloci` is not specified, all the loci present in the sequence
file specified by `seqfile` will be analysed.  
**DEPENDENCIES**  
`nloci` must be less than or equal to the number of loci available
in the sequence file specified by `seqfile`.  
**COMMENTS**  
If `nloci` is less than the number of loci present in the file
specified by `seqfile` only the first `nloci` present in that
file will be analysed.  
**EXAMPLES**
```
nloci = 1
nloci = 2000
```

### 18 model
------------------------------------------------------------------------
```
model = s [s]
```
**DESCRIPTION**  
Specifies the DNA substitution model (or amino acid substitution model)
that will be used. For DNA sequences the simplest (default) model that
can be chosen is the Jukes-Cantor (JC69) model and the most complex is
the General Time-Reversible (GTR) model. A range of models of
intermediate complexity are also available. For amino acid sequences,
most of the commonly used amino acid substitution matrix models are
available. If a model other than `Custom` is specified it will apply to
all loci. If `Custom filename` is specified, an additional file
`filename` must be provided that contains one or more lines with the
format:
```
loci_indices data_type model
```
where `loci_indices` is an index or range of indexes for loci (indexed
according to their order of occurrence in the sequence data file),
`data_type` is either `DNA` or `AA`, and model specifies the
substitution model variable for the indexed loci as described above. For
example, the entry
```
1-10 DNA JC69
11 AA DAYHOFF
```
specifies that loci 1 to 10 comprise DNA sequences that will have the
`JC69` model applied to them, while locus 11 is an amino acid sequence
that will have the `DAYHOFF` model applied to it.  
**VALUES**  
For DNA sequences: `JC69`, `K80`, `F81`, `HKY`,
`T92`, `TN93`, `F84`, `GTR`.  
For amino acid sequences: `DAYHOFF`, `LG`, `DCMUT`,
`JTT`, `MTREV`, `WAG`, `RTREV`, `CPREV`, `VT`,
`BLOSUM62`, `MTMAM`, `MTART`, `MTZOA`, `PMB`,
`HIVB`, `HIVW`, `JTTDCMUT`, `FLU`, and `STMTREV`.  
Among-locus model variation: `Custom filename`.  
**DEFAULT**  
`JC69`  
**DEPENDENCIES**  
If `model = GTR` the variables `Qrates` and `basefreqs` can
be optionally defined which provide priors on the rate matrix and base
frequencies, respectively (default priors will be used otherwise).  
**EXAMPLES**
```
model = GTR
model = CPREV
model = Custom models.txt
```

### 19 Qrates
------------------------------------------------------------------------
```
Qrates = b +f +f +f +f +f +f
```
**DESCRIPTION**  
Specifies whether the exchangeability parameters in the GTR model are
fixed or variable and either their fixed values or prior relative mean
values.  
**VALUES**  
`b`, specifies either fixed (1) or variable (0) exchangeability
parameters in the GTR model.  
`+f +f +f +f +f +f`, specify the exchangeability parameters in the
order $a, b, c, d, e, f$ for TC, TA, TG, CA, CG, and AG, as described in
[Yang 2014](https://global.oup.com/academic/product/molecular-evolution-9780199602605?lang=es&cc=gb).  
**DEFAULT**  
`0 1 2 1 1 2 1`  
**DEPENDENCIES**  
The variable `Qrates` is defined only if `model = GTR`.  
**COMMENTS**  
Two example cases are discussed here with either fixed, or variable
rates. The following specifies fixed rates:
```
Qrates = 1  2 1 1 1 1 2
```
This specifies fixed rates at $a = 2, b = c = d = e = 1$ and $f = 2$, so
that the transition rate is twice as high as the transversion rate, and
the model corresponds to K80 or HKY. Note that only the relative rates
matter, because the rate matrix ($Q$) is scaled so that the average rate
(over the base frequencies) is 1 and branch lengths are measured in the
expected number of mutations/substitutions per site. Nevertheless the
program expects six rate parameters.
```
Qrates = 0  10 5 5 5 5 10 
```
This specifies variable rates with a prior on $a, b, c, d, e$ that is a
Dirichlet distribution for every locus. The above specifies Dir(10, 5,
5, 5, 5, 10), with parameters
$\alpha_a = 10, \alpha_b = \alpha_c = \alpha_d = \alpha_e = 10$, and
$\alpha_f = 10$. The rates generated from the Dirichlet sum to 1, and
they are rescaled. Note that larger values for those parameters mean
less variance, while the mean for $a$, say, is given by
$\alpha_a/\alpha$ with
$\alpha = \alpha_a + \alpha_b + \alpha_c + \alpha_d + \alpha_e +
\alpha_f$. This specifies an average transition/transversion rate ratio
of 2 (if base frequencies are all fixed at $\frac{1}{4}$).  
**EXAMPLES**
```
Qrates = 1 2 1 1 1 1 2
Qrates = 0 1 1 1 1 1 1
```

### 20 basefreqs
------------------------------------------------------------------------
```
basefreqs = b +f +f +f +f
```
**DESCRIPTION**  
Specifies whether the base frequency parameters are fixed or variable
and either their fixed values or prior relative mean values.  
**VALUES**  
`b`, specifies either fixed (1) or variable (0) base frequency
parameters in the GTR model.  
`+f +f +f +f`, specify the base frequency parameters in the order
TCAG.  
**DEFAULT**  
`0 1 1 1 1`  
**DEPENDENCIES**  
The variable `basefreqs` is defined only if `model = GTR`.  
**COMMENTS**  
The parameters $\pi_T, \pi_C, \pi_A, \pi_G$ can be either fixed or
variable with a prior specified by the Dirichlet distribution with
$\alpha$ parameters specifying the relative mean frequencies.  
**EXAMPLES**
```
basefreqs = 1 0.15 0.35 0.15 0.35
basefreqs = 0 10 10 10 10
```

### 21 alphaprior
------------------------------------------------------------------------
```
alphaprior = f f +d
```
**DESCRIPTION**  
Specifies the prior probability density for the shape parameter
controlling the pattern of among-site rate variation.  
**VALUES**  
`f f +d`, specify the $\alpha_p$ and $\beta_p$ parameters of the
prior on the Gamma shape parameter and the number of rate categories
ncatG, respectively.  
**DEFAULT**  
$\alpha=\infty$ (no rate variation among sites)  
**COMMENTS**  
Among-site rate variation is assumed to follow a Gamma distribution with
a mean of 1 and shape parameter $\alpha$. The variance of rates among
sites is $1/\alpha$ and therefore a smaller $\alpha$ specifies more rate
variation. A discrete distribution is used to efficiently approximate
the Gamma distribution and ncatG is the number of rate categories used
in the approximation -- more rate categories provide a better
approximation but incur greater computational expense. A value for ncatG
of 4 is recommended. The prior on $\alpha$ is also a Gamma distribution
with parameters $\alpha_p$ and $\beta_p$. The variance and mean of the
prior on $\alpha$ are: 

$\textrm{Mean}(\alpha) = \frac{\alpha_p}{\beta_p}$,

$\textrm{Var}(\alpha) = \frac{\alpha_p}{\beta_p^2}$

**EXAMPLES**
```
alphaprior = 1 1 4
```

### 22 cleandata
------------------------------------------------------------------------
```
cleandata = b
```
**DESCRIPTION**  
Specifies whether columns in the alignment which have gaps or ambiguity
characters will be removed prior to analysis, or will instead be
retained for use in the likelihood calculation.  
**VALUES**  
`0`, specifies that columns in the alignment which have gaps or
ambiguity characters will be treated as missing data and used in the
likelihood calculation.  
`1`, specifies the removal of columns in the alignment which have
gaps or ambiguity characters.  
**DEFAULT**  
`0`  
**EXAMPLES**
```
cleandata = 1
cleandata = 0
```

### 23 thetaprior
------------------------------------------------------------------------
```
thetaprior = s [(f*, f f s)]
```
**DESCRIPTION**  
Specifies the distribution model and parameters of the prior
distribution on the contemporary and ancestral population parameters
$\theta$ and whether the MCMC integrates over the parameters or the
integration is instead performed analytically.  
**VALUES**  
`s`, specifies the form of the prior distribution on theta and
should be either `invgamma`, `gamma` or `beta` to specify either an
inverse gamma, gamma, or beta distribution, respectively.
`[(f*, f f s)]`, specifies the parameters of the prior distribution.
The form for the parameters depends on the type of distribution
specified in the first argument. The permissible combinations of
distributions and parameters are:
```
    invgamma a b
    invgamma a b e
```
which specifies an inverse gamma prior with $\alpha$ and $\beta$
parameters a and b, respectively. The first entry specifies that
$\theta$s are integrated over analytically and the second entry (with
trailing e) specifies that $\theta$s are estimated.
```
    gamma a b
```
which specifies a gamma prior with $\alpha$ and $\beta$ parameters a and
b, respectively.
```
    beta p q l u
```
which specifies a beta prior constrained to the interval $(l,u)$ with
parameters p and q, respectively. All the distributions above, except
the inverse Gamma distribution (without the e argument) integrate over
the $\theta$ parameters numerically as part of the MCMC, estimating the
posterior distribution of $\theta$ for each population.  
**COMMENTS**  
The inverse-gamma prior `invgamma a b`, where a and b are the parameters
$\alpha$ and $\beta$, has mean and variance: 

$\textrm{Mean}(\theta) = \frac{\beta}{\alpha - 1}$,

$\textrm{Var}(\theta) = \frac{\beta^2}{(\alpha - 1)^2 (\alpha - 2)}$.

For example, if $\alpha=3$ and $\beta=0.002$ the mean is
$0.002/(3 – 1) = 0.001$ (one variable site per kb on average). Note that
all $\theta$ parameters in the MSC model (for both modern species and
extinct ancestral species) are assigned the same specified prior
distribution with identical parameters. The inverse-gamma is a conjugate
prior for $\theta$ (Hey and Nielsen, 2007), which means that both the
prior and the posterior of $\theta$ will be inverse-gamma. Use of the
conjugate prior allows the $\theta$ parameters to be integrated out
analytically, and thus the dimension of the parameter space is reduced.
This typically leads to improved mixing of the MCMC. However, with this
option, the posterior of the $\theta$ parameters will not be produced.
To estimate the $\theta$ parameters when using an inverse gamma prior,
add the letter e (or E) on the line, as follows
```
    thetaprior = 3 0.002 e
```
Whether $\theta$s are integrated out or estimated, other results (such
as the posterior probability of species trees or species-delimitation
models) should be identical if the same prior is used. The analytical
integration of $\theta$ is only possible with an inverse gamma prior
distribution, if any other prior is used the $\theta$s are automatically
estimated and the E flag is not needed and should not be used. The gamma
prior `gamma a b`, where a and b are the parameters $\alpha$ and
$\beta$, has mean and variance: 

$\textrm{Mean}(\theta) = \frac{\alpha}{\beta}$,

$\textrm{Var}(\theta) = \frac{\alpha}{\beta^2}$.

For example, if $\alpha=0.001$ and $\beta=1$ the mean is
$0.001/1 = 0.001$ (one variable site per kb on average). The beta prior
`beta p q l u` has mean and variance: 

$\textrm{Mean}(\theta) = \frac{p u + q l}{p + q}$,

$\textrm{Var}(\theta) = \frac{p q (u - l)^2}{(p + q)^2 (p + q + 1)}$.

For example, with `beta 2 18 1e-6 0.1` the mean is: 

$\frac{2(0.1)+18(10^{-6})}{2+18} = 0.0100009$ 

which is 10 variable sites per kb on average.  
**EXAMPLES**
```
thetaprior = invgamma 3 0.002
thetaprior = invgamma 4 0.001 e
thetaprior = beta 2 18 1e-6 0.1
thetaprior = beta 2 18 1e-6 0.01
thetaprior = beta 1 10 0.0001 0.2
thetaprior = gamma 0.001 1
```

### 24 tauprior 
------------------------------------------------------------------------
```
tauprior = s f f
```
**DESCRIPTION**  
Specifies the prior probability distribution of the divergence time
parameter for the root in the species tree.
**VALUES**\
`s f f`, specifies either the gamma or inverse gamma prior
distribution and the $\alpha$ and $\beta$ parameters, respectively for
$\tau_0$, the divergence time parameter for the root in the species
tree.
```
    gamma a b
    invgamma a b
```
specifies a gamma distribution, or an inverse gamma distribution,
respectively with a and b to be the parameters $\alpha$ and $\beta$,
respectively.  
**COMMENTS**  
If the specified prior on $\tau_0$, the age of the root of the tree is
an inverse-gamma distribution IG($\alpha, \beta$) with $\alpha > 2$ the
mean and variance are: 

$\textrm{Mean}(\tau_0) = \frac{\beta}{\alpha - 1}$,

$\textrm{Var}(\tau_0) = \frac{\beta^2}{(\alpha - 1)^2 (\alpha - 2)}$.

If the prior is a gamma distribution the mean and variance are:

$\textrm{Mean}(\tau_0) = \frac{\alpha}{\beta}$,

$\textrm{Var}(\tau_0) = \frac{\alpha}{\beta^2}$.

The remaining species divergence times are generated from the uniform
Dirichlet distribution [@Yang2010 eq. 2]. As an example, suppose that
`tauprior = 3 0.03` specifies the prior for $\tau_0$, the divergence
time parameter for the root in the species tree. The mean of $\tau_0$ is
then $0.03/(3 – 1) = 0.015$ (which means 1.5% of sequence divergence
between the root of the species tree and the present time). If the
mutation rate is $10^{–9}$ mutations/site/year, this distance will
translate to an absolute prior mean divergence time of 15 MY.  
**EXAMPLES**
```
tauprior = invgamma 3 0.03
tauprior = gamma 0.01 1
```

### 25 phiprior
------------------------------------------------------------------------
```
phiprior = f f
```
**DESCRIPTION**  
Specifies the parameters of the prior distribution on $\varphi$, the
introgression probability under the MSC-I model.  
**VALUES**  
`f f`, specify the parameters $\alpha$ and $\beta$, respectively, of
the beta distribution prior $\textrm{B}(\alpha,\beta)$ for $\varphi$.  
**DEPENDENCIES**  
If the `species&tree` variable `t` specifies a tree with
introgression nodes then `phiprior` must be defined.  
**COMMENTS**  
The introgression probability $\varphi$ takes values on the interval
$(0,1)$ and the prior distribution for $\varphi$ is assumed to be beta
$\textrm{B}(\alpha, \beta)$ which has mean and variance:

$\textrm{Mean}(\varphi) = \frac{\alpha}{\alpha+\beta}$,

$\textrm{Var}(\varphi) = \frac{\alpha \beta}{(\alpha + \beta)^2 (\alpha + \beta + 1)}$.

For example, `phiprior = 1 1` specifies the beta prior
$\textrm{B}(\alpha,\beta)$ with $\alpha = 1$ and $\beta = 1$ for $\varphi$, which has a mean of $1/2$ and a variance of $1/12$.  
**EXAMPLES**
```
phiprior = 1 1
phiprior = 1.5 1.5
```

### 26 locusrate
------------------------------------------------------------------------
```
locusrate = d+[(f f f s, s)]
```
**DESCRIPTION**  
Specifies the model of substitution rate variation among loci, as well
as the parameters, and the prior on the parameters of the model.  
**VALUES**  
`d+[(f f f s, s)]`, specifies the model and parameters. The variable
d+ takes one of 3 possible values with different numbers of additional
variables depending on d+. The permissible combinations are:

| **d+** | **[args]**                                                          | **Description**                    |
|--------|-----------------------------------------------------------------------|------------------------------------|
| 0      |                                                                       | Loci have identical rates                |
| 1      | f f f s                                                               | Locus rates variable and estimated |
|        | $\alpha_{\bar{\mu}} \,\, \beta_{\bar{\mu}} \,\, \alpha_{\mu_i}$ prior |                                    |
| 2      | s                                                                     | Locus rates variable and specified |
|        | Filename of file containing locus rates                               |                                    |


where $\alpha_{\bar{\mu}}$ and $\beta_{\bar{\mu}}$ are the parameters of
the Gamma distribution prior on the average rate among loci $\bar{\mu}$
and $\alpha_{\mu_i}$ is the parameter of the prior on the rate at the
$i$th locus given $\bar{\mu}$. The prior variable must be either `iid`
for a conditional iid (hierarchical) prior on $\mu_i$ or `dir` for a
Gamma-Dirichlet prior. The prior variable is optional, if it is not
specified `iid` is used.  
**DEFAULT**  
`0`  
**DEPENDENCIES**  
If both the `prior` option for the variable `locusrate` and the `prior`
option for the variable `clock` are set, they should take the same
value.  
**COMMENTS**  
The setting `locusrate = 0` (default) means that all loci have the same
mutation rate. Setting `locusrate = 1 f f f s` specifies a model with
rate variation among loci in which rates are estimated from the sequence
data. Including the first integer (1) which species the model there are
5 arguments, the last one (prior) being optional. Parameters
$\alpha_{\bar{\mu}}$ and $\beta_{\bar{\mu}}$ specify the shape and rate
parameters of the gamma distribution for the mean rate across loci
($\bar\mu$), while $\alpha_{\mu_i}$ and prior are used to specify the
locus rates ($\mu_i$) given the mean rate ($\bar\mu$). The option prior
can take two values dir (for Gamma-Dirichlet rates for loci), and iid
(for conditional i.i.d. or hierarchical prior).
`locusrate = 2 LocusRateFileName` specifies the fixed-rates model of
locus-rate variation ([Burgess and Yang 2008](https://doi.org/10.1093/molbev/msn148)). This is the strategy used by
[Yang 2002](https://doi.org/10.1093/genetics/162.4.1811), with the relative rates estimated by the distance to an
outgroup species. The relative locus rates are listed in the file: there
should be as many numbers in the file, separately by spaces or line
returns, as the number of loci (nloci). The program re-scales those
rates so that the average across all loci is 1 and then use those
relative rates as fixed constants. Specifically the mean rate across
loci ($\bar{\mu}$) is assigned a gamma prior:

$\bar{\mu} = \mathrm{G}(\alpha_{\bar{\mu}},\beta_{\bar{\mu}}),$

with mean and variance:

$\textrm{Mean}(\bar{\mu}) = \frac{\alpha_{\bar{\mu}}}{\beta_{\bar{\mu}}}$,

$\textrm{Var}(\bar{\mu}) = \frac{\alpha_{\bar{\mu}}}{\beta_{\bar{\mu}}^2}$.

When there are no fossil calibrations in the species tree, the rates
should all be relative. In this case we suggest fixing
$\alpha_{\bar{\mu}}=\beta_{\bar{\mu}}=0$ (with 0 causing the program to
fix both parameters at $\infty$), and the program will fix the mean rate
across loci at $\bar\mu = 1$. Otherwise one can use equal and large
values such as $\alpha_{\bar{\mu}}=100$ and $\beta_{\bar{\mu}}=100$, so
that $\bar\mu$ is nearly fixed at 1. Given the mean rate $\bar\mu$, two
priors are available to specify the locus rates $\mu_i$, with
$i = 1, 2, \cdots, L$, where $L$ is the number of loci (`nloci`). If
`prior = dir` (for Gamma-Dirichlet distribution of locus rates), the
total rate $L\bar\mu$, given the mean rate ($\bar\mu$), is partitioned
into locus rates ($\mu_i$), using the concentration parameter
$\alpha_{\mu_i}$ The model and notation follow [Burgess and Yang 2008](https://doi.org/10.1093/molbev/msn148) [eq. 4] and
[Dos Reis et al 2014](https://doi.org/10.1093/sysbio/syu020) [eqs. 3-5]. If `prior = iid` (for conditional-i.i.d. or
hierarchical prior of locus rates), the locus rates ($\mu_i$) are
i.i.d. given the mean rate ($\bar\mu$):

$\mu_i \sim \mathrm{G}(\alpha_{\mu_i},\alpha_{\mu_i}/\bar{\mu})$. 

This model is described in [Zhu et al 2015](https://doi.org/10.1093/sysbio/syu109) [eq. 8] and is also implemented in
[mcmctree]{.smallcaps}. In both the Gamma-Dirichlet and the conditional
i.i.d. models, parameter $\alpha_{\mu_i}$ is inversely related to the
extent of rate variation among loci, with a large $\alpha_{\mu_i}$
meaning similar rates among loci. If all loci are noncoding, the rates
are probably similar, so $\alpha_{\mu_i} =$ 10 or 20 may be reasonable,
while for coding loci or exons, $\alpha_{\mu_i} =$ 2 or 1 may be
appropriate. The $\alpha_{\mu_i}$ parameter may affect the estimates of
the population size parameter ($\theta$) for the root node on the
species tree ([Burgess and Yang 2008](https://doi.org/10.1093/molbev/msn148)). The $L$ locus rates ($\mu_i$) are
parameters in the model. If $\alpha_{\bar{\mu}} > 0$, the mean rate
($\bar\mu$) is a parameter as well.  
**EXAMPLES**
```
locusrate = 0
locusrate = 2 rates.txt
locusrate = 1 0 0 2
locusrate = 1 2 3 2 dir
locusrate = 1 1 1 1 iid
```

### 27 clock
------------------------------------------------------------------------
```
clock = +d [f f f s s]
```
**DESCRIPTION**  
Specifies whether a strict molecular clock (equal substitution rates
among lineages) is used or instead a specific variable clock model
(allowing variation in substitution rates among lineages). If a variable
clock model is specified, the priors on the parameters of the variable
clock model are also specified.  
**VALUES**  
`+d`, specifies a strict clock (1), a variable clock with
independent rates among branches (2), or a variable clock with
autocorrelated rates between ancestral and descendent branches (3).  
`f f f s s`, are required when `clock =` 2 or 3 and specify
respectively the parameters $\alpha_{\bar{\nu}}$, $\beta_{\bar{\nu}}$
and $\alpha_{\nu_i}$, the prior distribution for the locus rate
($\mu_i$) which is either `dir` for Dirichlet or `iid` for conditional
i.i.d (hierarchical) prior, and the distribution for the branch rates
($r_{ij}$) given the locus rate ($\mu_i$) which is either `G` for Gamma
or `LN` for log-normal. Parameters $\alpha_{\bar{\nu}}$ and
$\beta_{\bar{\nu}}$ are the parameters of the gamma distribution for the
average of variances across loci ($\bar\nu$), while $\alpha_{\nu_i}$ and
`prior` are used to specify the variance ($\nu_i$) for locus $i$ given
the average ($\bar\nu$). The specification of the distribution of
$\nu_i$ given $\bar\nu$ follows the same procedure as the specification
of the distribution of $\mu_i$ given the mean rate $\bar\mu$ when
modeling among-locus rate variation; see notes above about the
`locusrate` variable.  
**DEFAULT**  
`1`  
**DEPENDENCIES**  
If both the `prior` option for the variable `locusrate` and the `prior`
option for the variable `clock` are set, they should take the same
value.  
**COMMENTS**  
As an example, the specification (`clock = 2 10.0 100.0 5.0 iid G`)
specifies a local-clock model in which the mutation rate drifts over
branches, independently among loci. A detailed explanation is as
follows. First $\bar\nu$ is assigned a gamma distribution
$G(10.0, 100.0)$, with mean 0.1. Given $\bar\nu$, the conditional
i.i.d. prior means that
$\nu_i \sim G(\alpha_{\nu_i}, \alpha_{\nu_i}/\bar\nu)$ with shape
parameter $\alpha_{\nu_i} = 5.0$ and mean $\bar\nu$, for
$i = 1, 2, \cdots, L$. Note that $\alpha_{\nu_i}$ is inversely related
to the variance of $\nu_i$: use small values of $\alpha_{\nu_i}$ (2, 1,
or 0.5) if you believe that $\nu_i$ varies among loci (meaning that the
clock nearly holds at some loci but is seriously violated at others).

As another example, the specification
(`clock = 2 10.0 100.0 5.0 dir LN`) means that $\bar\nu$ is assigned a
gamma prior $G(10.0, 100.0)$, and then the sum $L\bar\nu$ is partitioned
into $\nu_i$ (for $i = 1, 2, \cdots, L$) according to the Dirichlet
distribution (`prior = dir`) with concentration parameter
$\alpha_{\nu_i} = 5.0$. Again large $\alpha_{\nu_i}$ means the same
extent of clock violation at different loci.

Given the locus rate $\mu_i$ (specified using the `locusrate` variable)
and the variance $\nu_i$ (specified using the `clock` variable) for
locus $i$, the different lineages may have different rates at the locus,
and those rates are independent among loci. If distribution specified is
`G` (for gamma), the rate for species-tree branch $j$ at locus $i$ has
the following gamma distribution

$r_{ij} | \mu_i, \nu_i \sim \mathrm{G}(\mu_i^2/\nu_i, \mu_i/\nu_i)$.

This has mean $\mu_i$ and variance $\nu_i$. Alternatively if
distribution = `LN` (for log-normal), the rate for species-tree branch
$j$ at locus $i$ has the following the log-normal distribution

$r_{ij} | \mu_i, \nu_i \sim \mathrm{LN}(\mu_i, \nu_i)$,

where $\mu_i$
is the mean of the LN distribution and $\nu_i$ is the variance parameter
of the lognormal. Note that the rate-drift model specifies rates for
branches on the species tree (rather than on the gene tree) for each
locus, and gene-tree branches residing in the same population or species
have the same rate. For example, if all sequences at a locus are from
the same species and all coalescent events occur in that species (before
reaching an ancestral species), all branches on the gene tree will have
the same rate even if the relaxed-clock model allows rates among
species. In contrast, if a gene-tree branch passes several species or
populations, the different segments of the branch will have different
rates. The branch length on the gene tree is calculated by summing up
the lengths of those segments (with the length of each segment being the
product of the rate and the time duration for the segment).

The option `clock = 3` is similar to clock 2 but specifies the
autocorrelated-rates model. `clock = 3` with the log-normal distribution
(LN) specifies the geometric Brownian motion model of [Rannala and Yang 2007](https://doi.org/10.1080/10635150701420643). This
assigns a rate to each species-tree branch, that is, to the mid-point of
the branch. Given the rate at the species-tree root ($\mu_i$ at locus
$i$), the rates for the two branches around the root are specified. Then
given the rate for each ancestral branch, the rates for its two daughter
branches are specified, by integrating over the rate at the internal
node that is ancestral to the daughter branches. See figure 1 and
equations 3-8 in [Rannala and Yang 2007](https://doi.org/10.1080/10635150701420643). The rates for all species-tree branches
are thus assigned through a pre-order tree traversal, starting from the
root moving to the tips, until all branches are visited. If `clock = 3`
is specified with the gamma distribution (option `G`), the model works
as follows, using a similar pre-order tree traversal. First the two
branches at the species-tree root have the gamma distribution with mean
$\mu_i$ and variance $\nu_i$. Then given the rate for each ancestral
branch, the rates for its two daughter branches are specified as
independent gamma variables with the mean to be the rate of the parental
branch and with the variance to be $\nu_i$. Clock 3 is currently
implemented for the MSC model only and is unavailable under the MSC-I
model.

Note that a larger variance $\nu_i$ implies a more serious violation of
the molecular clock at locus $i$. Also note that $\nu_i$ will be similar
to $\bar\nu$, especially if $\alpha_{\nu_i}$ is large, and $\bar\nu$ has
prior mean $\alpha_{\bar{\nu}} / \beta_{\bar{\nu}}$. If
$\alpha_{\bar{\nu}} = 10$ and $\beta_{\bar{\nu}} = 100$, the prior mean
will be 0.1. Thus, for the log-normal model $\nu = 0.5$ represents a
serious violation of the clock while $\nu < 0.1$ represents only a
slight violation.  
**EXAMPLES**
```
clock = 1
clock = 2 10.0 100.0 5.0 dir LN
clock = 3 10.0 50.0 3.0 dir G
```

### 28 heredity
------------------------------------------------------------------------
```
heredity = +d [(f f, s, )]
```
**DESCRIPTION**  
Specifies whether proportional differences exist between all $\theta$
values for different loci, and if so whether the proportionality
multipliers (called inheritance scalars) are estimated from the data, or
set to specified fixed values. If inheritance scalars are estimated, the
parameters of the prior also must be specified.  
**VALUES**  
`+d`, specifies whether no differences in $\theta$ exist among loci
(0), proportional differences exist and inheritance scalars are
estimated (1), or proportional differences exist and inheritance scalars
are fixed to specified values (2).  
`f f`, is only specified when `+d` is 1 and specifies the $\alpha$
and $\beta$ parameters, respectively, of the Gamma prior on the scalar
$s_i$ for locus $i$.  
`s`, is only specified when `+d` is 2 and specifies the name of a
file containing $L$ (number of loci) heredity scalars (which must be
positive numbers) separated by spaces or linebreaks.  
**DEFAULT**  
`0`  
**COMMENTS**  
`heredity = 0` is the default and constrains $\theta$ to be equal across
loci. `heredity = 1` or `2` specifies two models that allow $\theta$ to
vary among loci, which may be useful for analyses combining data from
autosomal loci, loci on sex chromosomes, and/or mitochondrial loci. Such
mixed data, may be expected to have proportional differences in
effective population sizes among loci so that inheritance scalars
[Hey and Nielsen 2004](https://doi.org/10.1534/genetics.103.024182) should be applied. For example, in diploid sexual organisms
with strict maternal inheritance of mtDNA we expect the effective
population size of a mitochondrial locus to be $1/4$ that of a nuclear
locus. Other factors such as natural selection may also cause $\theta$
to deviate from the neutral expectation even among autosomal loci.
BPP implements two options for allowing such variations.
The first option (`heredity = 1`) specifies that locus-specific
inheritance scalars $s_i$ be estimated, using a Gamma prior with
parameters $\alpha$ and $\beta$ specified by the user. The prior mean
and variance of $s_i$ are:

$\textrm{Mean}(s_i) = \frac{\alpha}{\beta}$,

$\textrm{Var}(s_i) = \frac{\alpha}{\beta^2}$.

For example, `(heredity = 1 4 4)` specifies a Gamma prior
$\textrm{G}(4, 4)$, with mean $4/4 = 1$ and variance $1/4$, for the
inheritance scalar at each locus. The MCMC then generates the posterior
distribution of $s_i$ for each locus. The second option (`heredity = 2`)
allows the user to specify the $s_i$ in a file so they remain fixed
constants during the MCMC run.  
**EXAMPLES**
```
heredity = 0
heredity = 1 4 4
heredity = 2 scalars.txt
```

### 29 checkpoint
------------------------------------------------------------------------
```
checkpoint = +d [(+d)]
```
**DESCRIPTION**  
Specifies whether one or more checkpoint files are created to allow the
BPP program to resume using the current state of the MCMC at the
iteration during which the checkpoint file was created.  
**VALUES**  
`+d[(+d)]`, specifies the iteration at which the first checkpoint
file is created (first argument) and how frequently checkpoint files are
subsequently created (second optional argument).  
**DEFAULT**  
If the checkpoint value is undefined no checkpoint file is created.  
**COMMENTS**  
The checkpoint files are named OUTFILE.Z.chk where outfile is the name
specified by the `outfile` variable and Z is the number of the
checkpoint file (the first checkpoint file is labeled 1, the second 2,
and so on). To resume execution of the program starting from a
checkpoint use the BPP `–resume` switch followed by the checkpoint file
name. For example, `bpp –resume outrun1.1.chk`.  
**EXAMPLES**
```
checkpoint = 100000
checkpoint = 100000 10000
```

### 30 constraintfile
------------------------------------------------------------------------
```
constraintfile = s
```
**DESCRIPTION**  
Specifies the name of a file containing information for placing
topological constraints on inferred trees and specifying outgroups for
tree rooting.  
**VALUES**  
`s`, a string specifying the name of the constraint file.  
**DEFAULT**  
If `constraintfile` is undefined no constraints or outgroups are used in
the analysis.  
**EXAMPLES**
```
constraintfile = myconstraints.txt
```

### 31 threads
------------------------------------------------------------------------
```
threads = +d [+d +d]
```
**DESCRIPTION**  
Specifies the number of CPU threads to be used and optionally the
specific threads to be used.  
**VALUES**  
`+d [+d +d\]`, either 1 or 3 positive integers. A single integer
argument specifies the number of threads to use. If three integers are
given the first specifies the number of threads, the second specifies
the index of the first thread used and the third specifies the
increment. For example, `12 4 1` specifies that 12 threads should be
used, the first with index 4 and incrementing by one so that threads
4-15 are used.  
**DEFAULT**  
If `threads` is undefined, the calculations for each locus are placed on
different threads when possible.  
**DEPENDENCIES**  
The number of specified threads cannot exceed the number of loci
specified.  
**EXAMPLES**
```
threads = 4
threads = 8 4 1
```

## Example control file

To examine the structure of a typical BPP control file, we consider the
example file `A00.bpp.ctl` contained in the BPP distribution
subdirectory `examples/frogs`. The goal here is to provide an example of a
correctly formatted control file and briefly explain the syntax and
meaning of specified variables, as well as proposing some conventions
for writing optional comments in control files. An exhaustive
description of all control file variables, their effects, and
permissible values is provided in the definitions of
[BPP Control File Variables](#bpp-control-file-variables). 
In depth descriptions of specific model and prior choices can be found in
[Substitution Models](#substitution-models), [Introgression and Migration Models](#introgression-and-migration-models), and [Isolation with Migration Models](). 
Control files illustrating the 4 methods of analysis A00, A01, A10 and A11 are found in
[Methods of Analysis](#methods-of-analysis). The content of the control file `A00.bpp.ctl` is displayed below:
```
    seed =  -1

    seqfile = frogs.txt
    Imapfile = frogs.Imap.txt
    outfile = out.txt
    mcmcfile = mcmc.txt

    # fixed species tree and delimitation
    speciesdelimitation = 0 
    speciestree = 0

    species&tree = 4  K  C  L  H
                   9  7 14  2
                   (((K, C), L), H);

    # sequence data are unphased for all 4 populations
    phase =   1  1  1  1

    # 0: no data (prior); 1:seq like
    usedata = 1  

    # number of data sets in seqfile
    nloci = 5  

    # remove sites with ambiguity data (1:yes, 0:no)?
    cleandata = 0    

    # invgamma(a, b) for theta
    thetaprior = invgamma 3 0.004 E  

    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's
    tauprior = invgamma 3 0.002    

    *     heredity = 1 4 4
    *     locusrate = 1 5

    # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
    finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  

    # MCMC samples, locusrate, heredityscalars, Genetrees
    print = 1 0 0 0   

    burnin = 8000

    sampfreq = 2

    nsample = 100000
```
**A comment about comments**  
Optional comments may be added to the control file. If either a '\#' or
a '\*' symbol occurs on a line, everything following the symbol up to
the end of the line will be ignored by the BPP program. This allows
optional comments to be added by the user either to explain choices (or
meanings) of variables, or to retain optional alternative variable
settings for future use. The control files presented in this
documentation use the different comment-initiating symbols for distinct
purposes:

1.  '\#' initializes a comment that lists and/or explains the possible
    settings of the variable found immediately below the comment.

2.  '\*' initializes a comment that retains deactivated variable
    settings (for future use). If an alternative to a defined variable
    it follows that variable.

It is good practice to make use of these conventions for comments in
your own bpp control files.

### Control file variables explained

Here we walk through the variables in the above control file example
explaining the meaning of each and the valid values that can be
specified. An exhaustive description of all control file variables and
permissible values is given in [Control File](#control-file).
```
    seed = -1
```
This line specifies the seed value used to initialize the random number
generator. Because BPP uses pseudorandom numbers to drive a stochastic
algorithm for inference the results will vary between analyses initiated
with different seeds, but will be identical if the same seed is used
(because pseudorandom numbers are generated using a deterministic
algorithm started with the value of seed). Setting `seed = -1` specifies
that BPP use the clock on the computer to generate a seed. With this
setting, every run of the program will automatically be started with a
different seed. The seed that was used for a run is stored in a file
generated by the program called SeedUsed.
```
    seqfile = frogs.txt
    Imapfile = frogs.Imap.txt
    outfile = out.txt
    mcmcfile = mcmc.txt
```
These lines specify the names of the files used to run the program and
summarize the output. The program will assume that all these files are
located in the current working directory since no path is specified. The
line `seqfile = frogs.txt` specifies that the sequence data are in the
file named `frogs.txt`. See [Sequence File](#sequence-file) for a 
description of the sequence data file format.
The line `Imapfile = frogs.Imap.txt` specifies that the map file is
named `frogs.Imap.txt`. See [Imap File](#imap-file) for
a description of the map file format. The lines `outfile = out.txt` and
`mcmcfile = mcmc.txt` specify that the program output and mcmc samples
are printed to files named `out.txt` and `mcmc.txt`, respectively. The
format of the output file depends on the type of analysis being
performed and is described in [Methods of Analysis](#methods-of-analysis).
```
    # fixed species tree
    speciesdelimitation = 0 
    speciestree = 0
```
The line `speciesdelimitation = 0 ` specifies that the number of species
is fixed (not delimited) and the line `speciestree = 0` specifies that
the species phylogeny (topology) is fixed (not estimated). This
combination specifies an analysis of type A00, which is parameter
estimation under a MSC using multiple-loci and multiple-species data on
a fixed species tree ([Rannala and Yang 2003](https://doi.org/10.1093/genetics/164.4.1645); [Burgess and Yang 2008](https://doi.org/10.1093/molbev/msn148)). Analysis A00 is
described in detail in [A00: Demography and Divergence Time Estimation](#a00-demography-and-divergence-time-estimation).
```
    species&tree = 4  K  C  L  H
                   9  7 14  2
                   (((K, C), L), H);
```
The line `species&tree = 4 K C L H` specifies the number of
species/populations and their names. This example specifies 4 species in
the data, which are K, C, L and H. The next line ` 9 7 14 2` specifies
the maximum number of sequences at any locus for each species/population
(in the same order as the species listing on the previous line). This
example specifies a maximum of 9 sequences for K, 7 for C, 14 for L and
2 for H. The program interprets these numbers as follows. First, if the
number of sequences for a species is 2 or greater, $\theta$ for that
species will be a parameter estimated by the program. Second, the sum of
the numbers of sequences over all species specifies the maximum number
of sequences at a locus. However, specifying 2 sequences for a species
when there are actually 10 sequences for that species at some loci in
the data file has no harmful effects. Note that if there are $s$ species
in the species tree, the model will minimally contain the following
parameters: $s – 1$ species divergence times ($\tau$s) and $s – 1$
ancestral species $\theta$s. A $\theta$ parameter will also exist for
each contemporary species with at least two sequences present at some
locus. If the data are for a single species/population, the model will
contain one parameter only, $\theta$ for that species.

The line `(((K, C), L), H);` specifies the (fixed) species tree in
Newick format ending with a semicolon (;). This example specifies a tree
with no introgression (the MSC model). Specification of the MSC-I model
with introgression is different (see
[Introgression and Migration Models](#introgression-and-migration-models) for a description of the MSC-I model). The
Newick and extended Newick formats for specifying phylogenetic trees are
described in [Introgression and Migration Models](#introgression-and-migration-models) and a good overview is found at
<https://en.wikipedia.org/wiki/Newick_format>.
```
    # sequence data are unphased for all 4 populations
    phase =   1  1  1  1
```
The `phase` variable is specified using a boolean (0/1) indicator
variable for each species/population, with 0 (the default) specifying
that sequences from that species are fully phased haplotypes and 1
specifying they are unphased diploid sequences. If this line is omitted,
all sequences are assumed to be fully phased (that is, value 0 for all
species). The program does not allow some sequences from a species to be
phased while other sequences from that same species are unphased. If the
variable is 1, each sequence from that species is treated as an unphased
diploid sequence, with heterozygous sites represented using ambiguity
characters `YRMKSW`. The chararcters `N`, `-` and `?` are allowed in the
sequence file and are all treated as `?`, but other ambiguities (HBVD)
are not allowed. Each sequence is then expanded/resolved into two
sequences, and the program averages over all the possible phase
resolutions according to the approach of [Gronau et al 2011](https://doi.org/10.1038/ng.937). For example a
sequence with two heterozygous sites R and Y is resolved into two
possible haplotypes: (i) A\...T and G\...C, and (ii) A\...C and G\...T.
Note that this option forces a sequence with no heterozygote sites to be
phased/resolved into two identical sequences.

The MSC model averages over different phase resolutions of the
heterozygote sites in the likelihood calculation. This option may
increase the memory usage and CPU time considerably if there are many
sequences with many heterozygote sites at one locus. The simulation
program MCcoal has been modified to simulate diploid unphased data as
well. This allows one to use simulations to examine the difference
between the analysis of the full data and the unphased data. If the
phase flag for a species is 0, all sequences from that species are
treated as resolved haplotype sequences, and ambiguities are interpreted
in the usual way (for example, Y means one nucleotide that is either T
or C).
```
    # 0: no data (prior); 1:seq like
    usedata = 1  
```
The line `usedata = 1` specifies that the likelihood is used in the MCMC
run to generate the posterior distribution of parameters. Setting this
to 0 allows one to run the MCMC without sequence data to generate the
prior (mostly used for debugging).
```
    # number of data sets in seqfile
    nloci = 5  
```
The line `nloci = 5` specifies the number of loci (alignments) to be 5.
There may be more loci in the sequence data file than is specified here.
For example, if 200 loci exist in the sequence data file and you specify
`nloci = 2` in the control file, BPP will read the first
two loci only.
```
    # remove sites with ambiguity data (1:yes, 0:no)?
    cleandata = 0    
```
The line `cleandata = 0` specifies that columns in the alignment which
have gaps or ambiguity characters are retained and used in the
likelihood calculation. Setting `cleandata = 1` instead specifies that
such data are removed prior to analysis.
```
    # invgamma(a, b) for theta
    thetaprior = invgamma 3 0.004 E  
```
The line `thetaprior = 3 0.004 E` specifies the inverse-gamma prior
IG($\alpha, \beta$) for the $\theta$ parameters, with the mean to be
$\beta/(\alpha-1)$. In the example, the mean is $0.004/(3 – 1) = 0.002$
(two substitutions per kb on average) and the option E specifies that
the program analytically integrates over $\theta$s. Note that all
$\theta$ parameters in the MSC model (for both modern species and
extinct ancestral species) are assigned the inverse-gamma prior with the
same parameters.

**Analytical integration of $\theta$ and the inverse gamma prior**  
The inverse-gamma is a conjugate prior for $\theta$ ([Hey and Nielsen 2007](https://doi.org/10.1073/pnas.0611164104))
which means that both the prior and the posterior of $\theta$ will be
inverse-gamma. Use of the conjugate priors allows the $\theta$
parameters to be integrated out analytically, and thus the dimension of
the parameter space is reduced. This typically leads to improved mixing
of the MCMC. However, with this option, the posterior of the $\theta$
parameters will not be produced. To estimate the $\theta$ parameters,
add the letter e (or E) on the line, as follows
```
    thetaprior = invgamma 3 0.002 e
```
Whether $\theta$s are integrated out or estimated, other results (such
as the posterior of the $\theta$ parameters or the posterior probability
of species trees or species-delimitation models) should be identical. A
useful strategy may be to run the A01, A10, and A11 analyses without
estimating the $\theta$ parameters, and after the MAP model (the best
species tree or the best species delimitation model) is identified, run
the A00 analysis with the model fixed to estimate all parameters.

Some notes about the inverse-gamma distribution. Since
BPP3.4, both the $\theta$ and $\tau$ parameters are
assigned the inverse gamma priors rather than the gamma priors in
version 3.3 or earlier. One difference is that the gamma is light-tailed
while the inverse-gamma is heavy-tailed, so that the inverse-gamma may
be less influential than the gamma if your prior mean is much too small.
The inverse-gamma distribution IG($\alpha, \beta$) has mean
$m = \alpha/(\beta-1)$ if $\alpha > 1$ and variance
$s^2 = \beta^2/[(\alpha-1)^2(\alpha-2)]$ if $\alpha > 2$, while the
coefficient of variation is $s/m = \sqrt{1/(\alpha-2)}$. If little
information is available about the parameters, you can use $\alpha = 3$
for a diffuse prior and then adjust the $\beta$ so that the mean looks
reasonable. For example, for the human species,
$\theta_H \approx 0.0006$, which means that two random sequences from
the human population are different at  0.06% of sites, less than 1
difference per kb. A sensible diffuse prior is then
`thetaprior = 3 0.002`, with mean $0.002/(3 – 1) = 0.001$.
```
    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's
    tauprior = invgamma 3 0.002    
```
The line `tauprior = 3 0.002` specifies the inverse-gamma prior
IG($\alpha, \beta$) for $\tau_0$, the divergence time parameter for the
root in the species tree. Other divergence times are generated from the
uniform Dirichlet distribution ([Yang and Rannala 2010, eq.2](https://doi.org/10.1073/pnas.0913022107)). In the example, the
mean is $0.002/(3 – 1) = 0.001$ (which means 0.1% average sequence
divergence between the root of the species tree and the present time).
If the mutation rate is $10^{–9}$ mutations/site/year, this distance
will translate to a divergence time of 1 MY.
```
    # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
    finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  
```
The line `finetune = 1: 5 0.001 0.001 0.001 0.3 0.33 1.0 ` specifies the
step lengths used in proposals for the parameters of the MCMC algorithm.
These step lengths can affect both mixing and convergence of the MCMC.
The first value, before the colon (1 in the example), is a boolean
switch, with 1 specifying automatic adjustments of step lengths by the
program and 0 specifying no automatic adjustments. Following the colon
are the step lengths for the proposals used in the program. If you
choose to let the program adjust the step lengths, burnin has to be
\>200, and the step lengths specified here will then be the initial step
lengths, with the program adjusting/optimizing step lengths using the
information collected during the burnin step. Automatic adjustment works
well in most cases and is usually the best option. Some notes about
manually adjusting the step lengths are provided below in section ?.
```
    # MCMC samples, locusrate, heredityscalars, Genetrees
    print = 1 0 0 0   

    burnin = 8000

    sampfreq = 2

    nsample = 100000
```
The line `print = 1 0 0 0` specifies the information that is printed out
during the MCMC run. It has several boolean variables (0 for off and 1
for on) that control the printouts during the MCMC. The first variable
specifies whether the MCMC samples are written into the main sample file
(mcmc.txt). The second flag is for locus-specific rate parameters (rate
$\mu_i$, variance parameter $\nu_i$ and species-tree branch rates for
locus $r_{ij}$). The third flag is for locus-specific heredity scalars
if those are estimated from the data. The fourth flag is for printing
locus-specific gene trees. The fifth flag is for printing locus-specific
substitution-rate parameters (qmat for the $Q$ matrix, freqs for base
frequencies, and alpha for gamma rates for sites).

MCMC samples are taken after the burnin (8000 iterations), and in this
example, are taken every 2 iterations, with a total of 100,000 samples
taken. The total number of MCMC iterations is `burnin` + `sampfreq`
$\times$ `nsample`. The resulting file can be large. For analysis A00
(`speciesdelimitation = 0, speciestree = 0`), this file can be analysed
using the R package or tracer. For other analyses (A01, A10, and A11),
the dimension of the sample is changing and cannot be analysed using
such traceplot programs.

**Combining results from different runs**  
The variable setting `print = -1` specifies that BPP will
not perform an MCMC analysis and will instead read an existing MCMC
sample file and summarize the results. Thus, with setting `print = 1`,
the `mcmc.txt` file receives the MCMC program output, but with
`print = -1`, it instead provides the program input. When using this
option be careful not to overwrite files that you intended to keep. The
`print = -1` option is useful for combining MCMC samples from multiple
runs to produce the posterior summary. Suppose you run the same analysis
3 times in different directories -- you can merge the mcmc.txt files
from the three runs into one file and run the program using the
`print = -1` option to summarize the posterior for the combined sample.
Note that if you forced one or more of the jobs to terminate (using the
Unix kill command for example) the last line of each mcmc.txt file may
incomplete. If that is the case, you should delete the incomplete line.
The header lines of any files being concatenated should also be deleted.
Do not combine the MCMC samples from different types of analyses (i.e.,
using different data files and/or variable settings or priors).

## Substitution Models

BPP includes several models that allow more realistic representations of
the processes underlying the evolution of DNA and amino acid sequences.
Fundamental to all molecular phylogeographic analyses is the model of
the nucleotide or amino acid substitution process itself. DNA
substitution models can accomodate different rates of substitution
between nucleotides due to factors such as transition versus
transversion bias and unequal stationary nucleotide frequencies. BPP
allows the use of many different DNA base substitution and amino acid
substitution models. For DNA sequences, the models range from a simple
Jukes-Cantor model (computationally efficient) to a sophisticated
General Time-Reversible (GTR) model (computationally expensive but more
realistic). For closely related species the Jukes-Cantor model is often
adequate, while the GTR often performs better for more distantly related
species. The remaining DNA substitution models (K80, F81, HKY, T92,
TN93, and F84) are all special cases (submodels) of the GTR model and
are of intermediate complexity, usually having one or two free
parameters ([see Yang 2014](https://global.oup.com/academic/product/molecular-evolution-9780199602605?lang=es&cc=gb)). Many fixed substitution matrix models (no
free parameters) for modeling substitutions among amino acid sequences
are also available.

Another factor to consider is variation in overall (or average)
substitution rates. The overall rate of substitution can vary in
multiple ways and, based on empirical evidence, this type of variation
is one of the most important factors to accommodate during
phylogeographic inference. There are 3 distinct types of substitution
rate variation that can be accommodated using the models implemented in
BPP: rate variation among sites, rate variation among loci, and rate
variation among species (or populations). Substitution rate variation
among sites, or among loci, is a common pattern and is likely due to the
pervasive effects of negative (purifying) selection. Variation of
substitution rates among species (or populations) is more often observed
among distantly related species and may reflect the effects of
differential selection, evolution of DNA repair mechanisms, or other
factors.

This Chapter briefly describes the models of DNA substitution, and of
substitution rate variation, implemented in BPP, the parameters and
priors associated with these models, and how to specify the different
models in the BPP control file.

### DNA Substitution Models

#### Jukes-Cantor model

One of the earliest and simplest models of DNA substitution, the
Jukes-Cantor model, assumes that all possible changes between the 4
nucleotides have equal rates and thus the stationary frequencies of
bases are $1/4$ each. This model is the default in BPP and will be used
if no model is specified in the control file. The control file variable
`model` specifies the substitution model and the explicit specification
of the JC69 model is:
```
    model = JC69
```
The JC69 model appears exceedingly unrealistic. However, for closely
related sequences (with few multiple substitutions) it often performs as
well as more sophisticated models while incurring much less
computational expense. There are no parameters associated with the JC69
model (other than the substitution rate) so no additional control file
variables associated with the JC69 substitution model need to be
specified.

#### GTR model

The General Time-Reversible (GTR) model implemented in BPP allows
different rates of substitution between different nucleotide pairs, with
the constraint that the process is "reversible" so that, for example,
the rate of transition from A to T equals the rate of transition from T
to A and so on. The GTR model accommodates factors such as
transition-transversion bias and should be used when sequences are more
divergent (distantly related species). The control file specification of
the GTR model is:
```
    model = GTR
```
Note that if the GTR model is used, two other variables may be specified
in the control file: `Qrates` and `basefreqs`. The `Qrates` variable
specifies whether the exchangeability parameters of the GTR model are
fixed or variable as well as either their fixed values, or prior mean
values, respectively. The control file specification is:
```
    Qrates = B a b c d e f
```
The variable `B` is either 0 (variable exchangeability parameters) or 1
(fixed exchangeability parameters). The variables `a b c d e f` are
positive numbers specifying the exchangeability parameters for TC, TA,
TG, CA, CG and AG substitutions respectively (see [Yang 1994](https://doi.org/10.1007/BF00178256)). Note that
TC is the rate of either T-\>C or C-\>T substitutions. Note that only
the relative rates matter because the rate matrix is scaled to have an
overall average substitution rate of 1. Thus, for the model of fixed
rates if x is any positive number a b c d e f and xa xb xc xd xe xf are
equivalent specifications. An example with fixed rates is as follows:
```
    Qrates = 1  2 1 1 1 1 2
```
which specifies $a=2,b=c=d=e=1,f=2$ so that the transition rate is twice
as high as the transversion rate, corresponding to a K80 or HKY model.
In most cases, users will want to use the GTR model with exchangeability
parameters variable. Using this specification, the posterior
distribution of the parameters will be estimated from the data. An
example with variable rates is as follows:
```
    Qrates = 0  10 5 5 5 5 10
```
The absolute magnitude of the numbers *does matter* in this case as it
determines the variance of the prior probability distribution on
exchangeability parameters. The variance is inversely proportional to
the absolute magnitude of the variables. The prior on the free
parameters $a, b, c, d, e$ is a Dirichlet distribution with parameters

$\alpha_a = 10, \alpha_b = \alpha_c = \alpha_d = \alpha_e = 5, \texttt{and} \, \alpha_f = 10$

The rates generated from the Dirichlet sum to 1 and are rescaled. If we
define

$\alpha = \alpha_a + \alpha_b + \alpha_c + \alpha_d + \alpha_e + \alpha_f$

then the mean for $a$ is $\alpha_a/\alpha$ and the variance for $a$ is

$\textrm{Var}(a) = \frac{\frac{\alpha_a}{\alpha}\left(1-\frac{\alpha_a}{\alpha}\right)}{\alpha + 1}$.

The means and variances of the other parameters are the same
(substituting subscripts). Thus, in this example, the mean for $a$ is
$\alpha_a/\alpha = 10/40 = 1/4$ and the variance is
$(1/4)(3/4)/41 = 4.6 \times 10^{-3}$. The prior specified above has an
average transition-transversion rate ratio of 2 (if base frequencies are
all fixed at 1/4).

The other control file variable that must be specified when using the
GTR model is the `basefreqs` variable which specifies whether the
stationary nucleotide (base) frequencies are fixed, or variable, and
either their fixed values, or relative prior means, respectively. The
control file specification is:
```
    basefreqs = B T C A G
```
where B is either 0 (variable base frequencies) or 1 (fixed base
frequencies) and T C A and G are numbers specifying the frequencies of
bases T, C, A and G, respectively. An example with fixed base
frequencies is as follows:
```
    basefreqs = 1  0.1 0.2 0.2 0.5
```
In most cases, users will want to use the variable base frequencies
setting which specifies that the posterior distribution of base
frequencies is estimated from the data. An example with variable base
frequencies is as follows:
```
    0  10 10 10 10
```
Here the 4 variables T C A G specify the alpha parameters of a Dirichlet
prior distribution on the nucleotide frequencies. So, the prior mean
frequency of base A in this example is $10/40 = 1/4$ and the variance is
$(1/4)(3/4)/41 = 4.6 \times 10^{-3}$. As noted above when discussing the
Dirichlet prior for exchangeability parameters the absolute magnitude of
the variables determines the variance of the Dirichlet prior. If a GTR
model is chosen and `Qrates` and `basefreqs` are not specified in the
control file default values of `Qrates = 0 1 2 1 1 2 1` and
`basefreqs = 0 1 1 1 1` will be used.

#### Other DNA substitution models

A range of models are available that fall between the JC69 model and the
GTR model in terms of complexity and number of parameters. The
additional models available include: K80, F81, HKY, T92, TN93, and F84.
The priors on the parameters of these models are pre-specified by the
program.

#### Amino acid substitution models

There are many amino acid substitution models available. These models
completely specify the substitution matrix so there are no free
parameters. Many of the models are described in [Yang 2014](https://global.oup.com/academic/product/molecular-evolution-9780199602605?lang=es&cc=gb). The available models are: DAYHOFF, LG, DCMUT, JTT, MTREV, WAG, RTREV, CPREV,
VT, BLOSUM62, MTMAM, MTART, MTZOA, PMB, HIVB, HIV,W, JTTDCMUT, FLU, and
STMTREV.

### Custom among-locus model variation

If a single substitution model is specified, all loci are assigned that
model. The control file option `model = Custom filename` can be used to
assign different models among loci, where the models are specified in
the file `filename`. This file should contain one or more lines with the
format:

    loci_indices data_type model

where `loci_indices` is an index or range of indexes for loci (indexed
according to their order of occurrence in the sequence data file),
`data_type` is either `DNA` or `AA`, and model specifies the
substitution model variable for the indexed loci as described above. For
example, the entry

    1-10 DNA JC69
    11 AA DAYHOFF

specifies that loci 1 to 10 comprise DNA sequences that will have the
`JC69` model applied to them, while locus 11 is an amino acid sequence
that will have the `DAYHOFF` model applied to it.

Models of Substitution Rate Variation
-------------------------------------

### Among-site rate variation

Substitution rates vary among sites due to many factors, including
purifying selection in coding regions. BPP allows among-site rate
variation to be accomodated using the canonical Gamma distributed rate
variation model ([see Yang 1994](https://doi.org/10.1007/BF00178256)). In this model, the mean rate among
sites is constrained to be 1 so that the expected number of
substitutions on a branch is proportional to its length (in units of
expected substitutions). If the mean of the Gamma distribution is fixed,
the only free parameter is the shape parameter $\alpha$. The variance of
rates among sites is $1/\alpha$ and therefore a smaller $\alpha$
specifies more rate variation. Gamma distributed rate variation is
specified using the control file variable `alphaprior`. If this variable
is not set the default is $\alpha = \infty$ so that the variance is zero
(rates are constant among sites). The prior on $\alpha$ is also a Gamma
distribution with parameters $\alpha_p$ and $\beta_p$. The control file
specification is as follows:

    alphaprior = a b ncatG

where $a = \alpha_p$ and $b = \beta_p$ are the parameters of the Gamma
distributed prior on the parameter $\alpha$ and `ncatG` is the number of
rate categories used to approximate the Gamma distribution. A discrete
distribution is used to efficiently approximate the Gamma distribution
and ncatG is the number of rate categories used in the approximation --
more rate categories provide a better approximation but incur greater
computational expense. A value for ncatG of 4 is recommended. The
variance and mean of the prior on $\alpha$ are: 

$\begin{aligned}
  \textrm{Mean}(\alpha) & = & \frac{\alpha_p}{\beta_p}, \nonumber \\
  \textrm{Var}(\alpha) & = & \frac{\alpha_p}{\beta_p^2} \nonumber\end{aligned}$
  
As an example, the setting `alphaprior = 1 1 4` specifies a prior on the
shape parameter $\alpha$ (with 4 rate categories) that has a mean and
variance of 1, thus the expected variance of $\alpha$ is also 1. This is
a fairly diffuse (uninformative) prior.

### Among-locus rate variation

Empirical evidence suggests that substitution rates can vary
substantially among loci, probably reflecting the influences of
selection. BPP incorporates a hierarchical prior to model rate variation
among loci. The models implemented in BPP allow each locus to have a
different substitution rate. For $L$ loci (alignments), there are $L$
substitution rate parameters, where $\mu_i$ is the rate at the $i$th
locus.

#### Prior distribution of the mean rate among loci

The mean rate across loci is

$\bar{\mu} = \frac{1}{L} \sum_{i=1}^L \mu_i$.

The prior distribution for the mean rate across loci $\bar{\mu}$ is a Gamma distribution with
parameters $\alpha_{\bar{\mu}}$ and $\beta_{\bar{\mu}}$ and mean and
variance, 

$\begin{aligned}
  \textrm{Mean}(\bar{\mu}) & = & \frac{\alpha_{\bar{\mu}}}{\beta_{\bar{\mu}}}, \nonumber \\
  \textrm{Var}(\bar{\mu}) & = & \frac{\alpha_{\bar{\mu}}}{\beta_{\bar{\mu}}^2}. \nonumber\end{aligned}$
  
If there are no fossil calibrations in the species tree, the rates
should all be relative. In this case, to avoid overparameterization we
suggest fixing the mean substitution rate across loci to be
$\bar\mu = 1$. This is done by setting
$\alpha_{\bar{\mu}} = \beta_{\bar{\mu}} = \infty$ (see discussion of
control file settings below), or setting both these parameters to large
values such as $\alpha_{\bar{\mu}} = \beta_{\bar{\mu}} = 100$, in which
case $\bar{\mu}$ is nearly fixed.

#### Prior distribution for locus-specific rates given mean rate among loci

Conditioning on the mean rate $\bar\mu$, two priors are available to
specify the probability distribution of locus-specific rates, $\mu_i$.
The control file variable argument `prior` (see discussion of control
file settings below) can take two values `dir` (for Gamma-Dirichlet
rates for loci), and `iid` (for conditional i.i.d. or hierarchical
prior). If `prior = dir` (Gamma-Dirichlet distribution of locus rates),
the total rate $L\bar\mu$, given the mean rate ($\bar\mu$), is
partitioned into locus rate variables , $\mu_i$, that are constrained to
sum to $L\bar\mu$ and have a Dirichlet distribution with concentration
parameter $\alpha_\mu$. The model and notation follow 
[Burgess and Yang 2008](https://doi.org/10.1093/molbev/msn148) eq.4
and [Dos Reis et al 2014](https://doi.org/10.1093/sysbio/syu020) eqs.3-5. 
If `prior = iid` (Conditional-i.i.d. or hierarchical prior of locus rates), the locus
rates, $\mu_i$, are independent and identically distributed i.i.d. given
the mean rate $\bar\mu$, so the prior distribution of $\mu_i$ is:

$\mathtt{P}(\mu_i) = \textrm{Gamma}(\alpha_{\mu_i}, \alpha_{\mu_i/\bar{\mu}})$.

This model is described in [Zhu et al 2015](https://doi.org/10.1093/sysbio/syu109) eq.8 and implemented in
[mcmctree]{.smallcaps}. In both the Gamma-Dirichlet and the Conditional
i.i.d. models, the parameter $\alpha_{\mu_i}$ is inversely related to
the extent of rate variation among loci, with a large $\alpha_{\mu_i}$
meaning similar rates among loci (less rate variation). If all loci are
noncoding, the rates are probably similar, so $\alpha_{\mu_I} =$ 10 or
20 may be reasonable, while for coding loci or exons,
$\alpha_{\mu_i} = 2$ or 1 may be appropriate. The $\alpha_{\mu_i}$
parameter may affect the estimates of the population size parameter
($\theta$) for the root node on the species tree ([Burgess and Yang 2008](https://doi.org/10.1093/molbev/msn148)). The $L$
locus rates ($\mu_i$) are parameters in the model. If
$\alpha_{\bar{\mu}} > 0$, the mean rate, $\bar\mu$, is a parameter as
well.

#### Control file settings for locusrate

The `locusrate` variable takes from 1 and 5 arguments. In the case of 5
arguments the form is shown below, with the last argument (prior) being
optional.

    locusrate = rate_variation a_mubar b_mubar a_mui prior 

The rate\_variation variable is either 0 (no rate variation, the
default), 1 (rate variation with rates inferred) or 2 (rate variation
with fixed user-specified rates). If rate\_variation = 0 no other
variables are specified, if rate\_variation = 2 then a second argument
specifies the name of a file containing $L$ specified fixed rates, one
for each locus. If rate\_variation = 1 the additonal variables `a_mubar`
and `b_mubar` specify the shape and rate parameters
($\alpha_{\bar{\mu}}$ and $\beta_{\bar{\mu}}$) of the gamma distribution
that is the prior for the mean rate across loci (`mubar` or $\bar\mu$),
and variable `a_mui` ($\alpha_{\mu_i}$) and `prior` are used to specify
parameter of the prior on the locus rates ($\mu_i$) given the mean rate
($\bar\mu$) and the form of that prior (`dir` or `iid`). Some example
configuration file settings are:

    #  (0: No rate variation among loci.  This is default)
    locusrate = 0
    # (1: estimate locus rates mui)
    locusrate = 1 10 10 5.0 iid
    # (1: estimate locus rates mui)
    locusrate = 1  0  0 5.0 iid
    # (2: locus rates from file)
    locusrate = 2 LocusRateFileName    

Note that `locusrate = 2 LocusRateFileName` specifies the fixed-rates
model of locus-rate variation ([Burgess and Yang 2008](https://doi.org/10.1093/molbev/msn148)). This is the strategy used
by [Yang 2002](https://doi.org/10.1093/genetics/162.4.1811), with the relative rates estimated by the distance to an
outgroup species. The relative locus rates are listed in the file: there
should be as many numbers in the file, separately by spaces or line
returns, as the number of loci (`nloci`). The program re-scales those
rates so that the average across all loci is 1 and then use those
relative rates as fixed constants. The following example nearly fixes
(first line), or fixes (second line) the mean rate $\bar{\mu}$. This is
the recommended practice when no fossil calibrations are available:

    # locusrate = 1 a_mubar b_mubar a_mui <prior>
    # (1: estimate locus rates mui)
    locusrate = 1 10 10 5 iid
    # (1: estimate locus rates mui, mubar=1 fixed)
    locusrate = 1  0  0 5 iid      

Note that `a_mubar = 0` and `b_mubar = 0` (with 0 meaning $\infty$)
cause the program to fix $\bar{\mu}$, while equal and large values such
as `a_mubar = 100` and `b_mubar = 100` cause cause $\bar\mu$ to be
nearly fixed at 1.

**The locus-rate prior in IMa.** The model of variable rates among loci
implemented here has some differences from a similar model implemented
in the IMa program [Hey and Nielsen 2004](https://doi.org/10.1534/genetics.103.024182). The biggest difference
appears to be the parameterization. BPP defines mutation
rate on a per-nucleotide basis, so the prior specifies that the
expectation of the mutation rate per site is constant among loci.
IMa defines mutation rate on a per-locus basis, so its
prior specifies that the expectation of the mutation rate per locus is
the same among loci. If locus one has 100 sites and locus two has 1000
sites, then IMa assumes that the per-site rate for locus
one is 10 times for locus two, while BPP assumes the same
per-site rate. Also IMa constrains the geometric mean of
rates across loci to be one, while BPP constrains their
arithmetic mean to be one.

### Among-species rate variation

The strict molecular clock model of evolution assumes that rates of
substitution for a particular gene, or genomic region, are the same in
different species. Early applications of relative rates tests suggested
that this assumption is often violated in empirical data, probably due
to factors such as selection and in more distantly related species
evolving DNA repair and proof-reading mechanisms.

One can allow rates of substitution to vary in an arbitrary way be
inferring unrooted phylogenetic trees but this does not allow species
divergence times to be estimated. Relaxed clock models allow the rate of
substitution to vary among contemporary and ancestral species but with
some constraints which allow species divergence times to be estimated.
The BPP program implements two relaxed clock models, as well as a strict
molecular clock model. The two relaxed clock models differ in terms of
their assumptions about rates in ancestor and descendent species. The
autocorrelated rates model assumes that descendents have rates that are
correlated with the rate in the ancestor, while independent rates model
assumes that the rates on each branch are iid from a common
distribution. For the independents rates model, the variance of the
distribution of rates among branches determines the degree of departure
from a molecular clock - a low variance means that the rates are almost
clocklike (vary little among branches).

#### Variation in the degree of departure from a molecular clock among loci

For a given set of species some loci may have rates that are highly
clock-like while other loci may have much greater variation in
substitution rates among species (branches of the species tree). To
accomodate this, the variance of the substitution rate among species is
allowed to vary among loci. The prior distribution of the average
variance $\bar{\nu}$ for all loci is

$P(\bar{\nu}) = \textrm{Gamma}(\alpha_{\bar{\nu}},\beta_{\bar{\nu}})$,

so that the expected average variance is
$\alpha_{\bar{\nu}}/\beta_{\bar{\nu}}$. If this value is larger loci
depart more from the molecular clock on average. Conditional on the
average variance parameter $\bar{\nu}$ the prior for the variance
$\nu_i$ at locus $i$ can be either a Gamma-Dirichlet distribution or a
Conditional i.i.d distribution (controlfile variables `dir` and `iid`,
respectively). Both distributions have a parameter $\alpha_{\nu_i}$ that
determines the variance of the variance in rates across loci. The
specification of the distribution of $\nu_i$ given the mean $\bar\nu$
follows the same procedure as the specification of the distribution of
the locus-specific substitution rate $\mu_i$ given the mean rate
$\bar\mu$ in the previous section; see notes above about the locusrate
variable. Note that $\alpha$ (control file variable `a_vi`) is inversely
related to the variance in $\nu_i$: use small values of $\alpha$ (2, 1,
or 0.5) if you believe that $\nu_i$ varies among loci (meaning that the
clock nearly holds at some loci but is seriously violated at others).

#### Relaxed clock with i.i.d rates among branches

The independent rates model of the relaxed molecular clock assumes that
for locus $i$ each branch has a rate that is independent and identically
distributed (i.i.d). The distribution has mean $\mu_i$ specified by the
`locusrate` model and variance $\nu_i$ specified by the prior
distribution for rates among loci discussed above. Two distributions are
available: a Gamma distribution or a log-normal distribution (control
file variables `LN` and `G`, respectively). Given the locus rate $\mu_i$
(specified using the `locusrate` variable) and the variance $\nu_i$
(specified using the `clock` variable) for locus $i$, the different
lineages may have different rates at the locus, and those rates are
independent among loci. If distribution = G (for gamma), the rate for
species-tree branch $j$ at locus $i$ has the following gamma
distribution

$r_{ij} | \mu_i, \nu_i \sim G(\mu_i^2/\nu_i, \mu_i/\nu_i)$.

This has mean $\mu_i$ and variance $\nu_i$.

Alternatively if distribution = `LN` (for log-normal), the rate for
species-tree branch $j$ at locus $i$ has the following the log-normal
distribution $$r_{ij} | \mu_i, \nu_i \sim \mathrm{LN}(\mu_i, \nu_i),$$
where $\mu_i$ is the mean of the LN distribution and $\nu_i$ is the
variance parameter of the lognormal.

#### Relaxed clock with autocorrelated rates among branches

The autocorrelated-rates model (control file option `clock = 3`) with
the log-normal distribution (LN) specifies the geometric Brownian motion
model of [Rannala and Yang 2007](https://doi.org/10.1080/10635150701420643). 
This assigns a rate to each species-tree branch,
that is, to the mid-point of the branch. Given the rate at the
species-tree root ($\mu_i$ at locus $i$), the rates for the two branches
around the root are specified. Then given the rate for each ancestral
branch, the rates for its two daughter branches are specified, by
integrating over the rate at the internal node that is ancestral to the
daughter branches. See figure 1 and equations 3-8 in [Rannala and Yang 2007](https://doi.org/10.1080/10635150701420643). The
rates for all species-tree branches are thus assigned through a
pre-order tree traversal, starting from the root moving to the tips,
until all branches are visited.

If the gamma distribution (G) is instead used the model works as
follows, using a similar pre-order tree traversal. First the two
branches at the species-tree root have the gamma distribution with mean
$\mu_i$ (`mu_i`) and variance $\nu_i$ (`nu_i`). Then given the rate for
each ancestral branch, the rates for its two daughter branches are
specified as independent gamma variables with the mean to be the rate of
the parental branch and with the variance to be $\nu_i$. The
autocorrelated relaxed clock model is only implemented for the MSC model
and is not available when using the MSC-I model.

**Rates vary among species tree branches**\
The rate-drift model specifies rates for branches on the species tree
(rather than on the gene tree) for each locus, and gene-tree branches
residing in the same population or species have the same rate. For
example, if all sequences at a locus are from the same species and all
coalescent events occur in that species (before reaching an ancestral
species), all branches on the gene tree will have the same rate even if
the relaxed-clock model allows rates to vary among species. In contrast,
if a gene-tree branch passes several species or populations, the
different segments of the branch will have different rates. In that
case, the branch length on the gene tree is the sum of the lengths of
those segments (with the length of each segment being the product of the
rate and the time duration for the segment).

#### Control file settings for clock

The `clock` variable takes 6 arguments when using a relaxed clock as
follows:
```
clock = clock_type a_vbar b_vbar a_vi prior dist
```
the argument `clock_type` takes on of three values: 1 (strict clock,
default), 2 (relaxed clock with i.i.d rates) and 3 (relaxed clock with
autocorrelated rates). In the case of the strict clock there are no
additional arguments. The next two variables `a_vbar` and `b_vbar`
specify the parameters ($\alpha_{\bar{\nu}}$ and $\beta_{\bar{\nu}}$) of
the Gamma prior specifying the average rate variance across loci (see
above). The variable `a_vi` specifies the variance parameter of the
prior on the locus specific variance, $\alpha_{\nu_i}$ (see above). The
variable `prior` is either `iid` or `dir` and specifies the prior on
$\nu_i$ given $\bar{\nu}$ (see above). The variable `dist` is either
`LN` or `G` specifying either a Gamma distribution of a log-normal
distribution for the branch rates conditional on the locus rate $\mu_i$
and among-branch variance $\nu_i$.

Several example control file entries for the `clock` variables are
provided below:
```
    # (1: strict clock, default)
    clock = 1 
    # (2: independent-rates)
    clock = 2 10.0 100.0 5.0 iid G
    # (3: correlated-rates)
    clock = 3 10.0 100.0 5.0 iid G
```
The specification (`clock = 2 10.0 100.0 5.0 iid G`) means the
following. First $\bar\nu$ is assigned a gamma distribution
$G(10.0, 100.0)$, with mean 0.1. Given $\bar\nu$, the conditional
i.i.d. prior means that $\nu_i \sim G(\alpha, \alpha/\bar\nu)$ with
shape parameter $\alpha = 5.0$ (`a_vi`) and mean $\bar\nu$, for
$i = 1, 2, \cdots, L$. As another example, the specification
(`clock = 2 10.0 100.0 5.0 dir LN`) means that $\bar\nu$ is assigned a
gamma prior $G(10.0, 100.0)$, and then the sum $L\bar\nu$ is partitioned
into $\nu_i$ (for $i = 1, 2, \cdots, L$) according to the Dirichlet
distribution (`prior = dir`) with concentration parameter $\alpha = 5.0$
(`a_vi`). Again large $\alpha$ means the same extent of clock violation
at different loci.

**Choosing $\alpha_{\bar{\nu}}$ and $\beta_{\bar{\nu}}$**  
A larger variance $\nu_i$ represents a greater violation of the
molecular clock at locus $i$. Note that $\nu_i$ will be similar to
$\bar\nu$, especially if $\alpha_{\bar{\nu}}$ (`a_vi`) is large, and
$\bar\nu$ has prior mean $\alpha_{\bar{\nu}}/\beta_{\bar{\nu}}$. If
(`a_vbar = 10 b_vbar = 100`), the prior mean will be 0.1. For the
log-normal model $\nu = 0.5$ suggests a serious violation of the clock
while $\nu < 0.1$ represents only a slight violation.

## Introgression and Migration Models

Another important factor influencing the results of species tree
inference is introgression (gene flow) between species (or populations).
BPP also implements a model of episodic introgression between species or
populations [Flouri et al 2020](https://doi.org/10.1093/molbev/msz296) as well as a model of ongoing gene flow (MSC-M).
Since version 4.1, BPP implements the MSC-I model
[Flouri et al 2020](https://doi.org/10.1093/molbev/msz296) and since version 4.4 it implements the MSC-M model.
Currently, introgression (MSC-I) and ongoing gene flow (MSC-M) analyses can
only be performed using the A00 method of analysis (fixed species tree
and delimitation, see [Methods of Analysis](#methods-of-analysis)). In other words, the user has to specify
the number of introgression events (or bands of gene flow), their
directions, and the populations involved on a fixed species tree. The
program will then estimate the parameters of the MSC-I (or MSC-M) models
using MCMC. We hope to implement MCMC moves that change the MSC-I (or MSC-M)
models and species tree topology in the future.

### The MSC-I Model

The MSC-I model involves three types of parameters: the species
divergence or hybridization times ($\tau$s), the population size
parameters ($\theta$s), and the introgression probabilities
($\varphi$s). The control file must specify an extended Newick format
tree with one or more introgression events as well as a prior
distribution for the introgression probability $\varphi$. The MSC-I model
is specified in the `species&tree` block using the extended Newick
notation ([Cardona et al 2008](https://doi.org/10.1186/1471-2105-9-532))
```
    species&tree = 3  A  B  C
    10 10 10
    ((A, (C)H[&phi=0.5,&tau-parent=yes])S, (H[&tau-parent=yes], B) T)R;
```
The rules for the extended Newick format are as follows. '(A, B)S'
specifies two branches from S to A and from S to B, while '(A)H'
specifies one branch from H to A. Every branch in the tree model is
represented once. Each tip species occurs once. Internal nodes for
speciation nodes may and may not be labeled, but hybridization nodes
must be labeled. The extended Newick notation is not unique.

Each hybrid node has two parental species. Models A-C are distinguished
using a meta-data variable `tau-parent`. The value yes means that the
parental population exists and has a separate age parameter $\tau$, with
an associated $\varphi$ parameter, while the value no means that such
parameters do not exist. For example, in model (A), $\tau_S$, $\tau_T$,
$\theta_{Hl}$ (for H-left), and $\theta_{Hr}$ (for H-right) are all
parameters; in model (B), $\tau_S$ and $\theta_{Hl}$ are not parameters
in the model while $\tau_T$ and $\theta_{Hr}$ are parameters; and
finally in model (C), none of those four parameters exists in the model.
Model (D) specifies a bidirectional introgression event between species
A and B.

![MSCI-models](https://github.com/bpp/bpp-manual/blob/191d107413566b145360597fdfdfef3e68ce9173/figures/fig-msci-models.png?raw=true)

**Four different types of MSC-I models implemented in BPP**  
The models are specified using the extended Newick format, as follows:
```
(A): ((A, (C)H[&phi=0.5,&tau-parent=yes])S, (H[&tau-parent=yes], B) T)R; 
(B): ((A, (C)H[&phi=0.5,&tau-parent=no])S, (H[&tau-parent=yes], B) T)R; 
(C): ((A, (C)H[&phi=0.5,&tau-parent=no])S, (H[&tau-parent=no], B) T)R; 
(D): ((A, (B) Y[&phi=0.3])X, (X[&phi=0.1])Y)R; 
(D): ((A, (B)Y)X, (X)Y)R; 
(D): ((A, Y)X, (B, X)Y) R; 
(D): ((A, y)x, (B, x) y) r;
```
The control file variable `phiprior` specifies the prior probability for
the $\varphi$ parameter, which is a beta distribution with parameters
$a$ and $b$. The syntax is
```
    phiprior = a b
```
where $a$ and $b$ are positive numbers. For example, `phiprior = 1 1`
specifies the beta prior beta($a, b$) with $a = 1$ and $b = 1$ for
$\varphi$, the introgression probability under the MSC-I model.

#### Autogeneration of extended newick graphs
Starting with BPP release 4.4.0, we have changed the definition of introgression probability ($\varphi$) so that it is assigned to the horizontal (introgression) branch.  This note illustrates the `msci-create` feature in bpp, which generates extended Newick notation for the MSC-I model from a data file which contains a binary species tree with introgression events specified using source and target branches on the tree.
The `--msci-create` option has been available since [BPP 4.2.1](https://github.com/bpp/bpp/issues/114).  The user prepares a input file (named `msci.txt` in our example), which includes a binary species tree in Newick format and specifies the introgression events by identifying the source and target branches involved in the introgression.  Run bpp using the following command:
```
$ bpp --msci-create msci.txt
```
This generates extended Newick notation for the MSC-I model, which can be copied into a bpp control file. 
The input file `msci.txt` uses four commands: `tree`, `define`, `hybridization`, and `bidirection`.

- `tree` defines the binary species tree, in Newick format.
- `define` defines an ancestral species or internal node label, as the most recent common ancestor of the tip species.  Alternatively you can label the internal nodes on the Newick tree.
- `hybridization` defines a hybridization/introgression event by specifying the source and target branches (which represent the source and target populations involved).
- `bidirection` defines a bidirectional introgression (BDI) event.

One may think of the binary species tree as describing the history of species divergences and add introgression events onto it as new (horizontal) branches.  The introgression probability is assigned to the newly created introgression branch. 

#### Examples using --msci-create
A series of examples follow. The introgression graph is displayed in a figure followed by the `msci=create` encoding and the resulting extended Newick format.
![Fig-MSCI2](https://github.com/bpp/bpp-manual/blob/36012225be26f3b82593c9a35eaaa3e7217a5fd0/figures/msci_example2.png?raw=true)

```
# Model A, version 1
tree (A,(B,C));
define T as B,C
define R as A,B
hybridization R A, T C as S H tau=yes, yes phi=0.10

# Model A, version 2
tree (A,(B,C)T)R;
hybridization R A, T C as S H tau=yes, yes phi=0.10

# The generated Newick notation for model A is 
# ((H[&phi=0.1,tau-parent=yes],A)S, (B,(C)H[&phi=0.9,tau-parent=yes])T)R;

# Model B1
tree (A,(B,C)T)R;
hybridization R A, T C as S H tau=no, yes phi=0.10

# The generated Newick notation for model B1 is 
# ((H[&phi=0.1,tau-parent=no],A)S, (B,(C)H[&phi=0.9,tau-parent=yes])T)R;

# Model B2
tree ((A,C)S,B)R;
hybridization R B, S C as T H tau=no, yes phi=0.10

# The generated Newick notation for model B2 is 
# ((A,(C)H[&phi=0.9,tau-parent=yes])S, (H[&phi=0.1,tau-parent=no],B)T)R;

# Model C
tree (A,(B,C)T)R;
hybridization R A, T C as S H tau=no, no phi=0.40

# The generated the Newick notation for model C is 
# ((H[&phi=0.4,tau-parent=no],A)S, (B,(C)H[&phi=0.6,tau-parent=no])T)R;

# Model D
tree (A,B)R;
bidirection  A R, B R as X Y phi=0.1,0.2
```
The above figure illustrates the specifications of models A, B, C, and D of [Flouri et al 2020](https://doi.org/10.1093/molbev/msz296).  In models A, B1, and C, the introgression event is from the source branch RA to the target branch TC, with a new introgression branch SH created.  The introgression probability ($\varphi = 0.1$) is assigned to the newly added introgression branch (SH) while $1 – \varphi$ is assigned to the other parent branch (TH).  Models A, B1, and C are distinguished by using the keyword `tau` (which means the same as `tau-parent` in the Newick notation of the MSci model). Model A assumes that the two parental species S and T of the hybridization node H have distinct ages from H, with $\tau_S > \tau_H$ and $\tau_T > \tau_H$.  Thus we have `hybridization R A, T C as S H tau=yes, yes phi=0.10`: the first `yes` means that parental node S on the source branch RA has a distinct tau from node H and that the second `yes` means that parental node T on the target branch TC has a distinct tau from node H.  Note that there are often multiple ways of specifying the same MSC-I model and here you can specify model A by starting with the binary tree `((A, C), B)` and adding branch TH as a introgression event.  Model A might appear to be nonsensical biologically since only contemporary species can exchange migrants.  Nevertheless, model A may be used to represent introgressions from a ghost species not sampled in the sequence dataset.  Suppose at time $\tau_S$ a speciation event generated two species SA and SH, and at $\tau_H$, species SH contributed migrants into species THC, but species SH has since become extinct or is otherwise not sampled in the data.  This scenario matches model A.
Model B1 assumes $\tau_S = \tau_H$ and $\tau_T > \tau_H$ so there is no new tau for parent S of the hybridization node H, and there is a new tau for parent T.  Thus we have `tau=no, yes`.  Model B2 works similarly.  In model C, none of parents S and T has a distinct age from H, so we have `tau=no, no`.  One interpretation is that species C is a hybrid species. In model D for the bidirectional introgression (BDI) model, there is no distinction of source and target branches.  We specify two phi values, assigned to the introgression (horizontal) branches: $\varphi_X = 0.1$ for node X (or into node X) and $\varphi_Y = 0.2$ for node Y (see above figure).

### The MSC-M Model
The multispecies-coalescent-with-migration model (or isolation-with-migration or IM model) is activated in BPP using the keyword migration.
A binary species tree is specified using an extended Newick format that includes labels for internal nodes; these are subsequently used in the control file to specify the source and target populations involved in migration.

#### Control file specification of MSC-M
The basic model is specified in the control file as follows, using a species tree for four species (A, B, C, D) as an example:
```
      ((A, B)S, (C, D)T)R;

     migprior = 2 200 
    migration = 2
                A C
                S C
```
Here S is the AB common ancestor, T is the CD common ancestor, while R is the ABCD ancestor.  Not all internal nodes need to be labeled but those involved in migration have to be.  The migration line specifies 2 migration connections: one connecting source population A to target population C, and another connecting source population S to target population C, with migration rates $M_{AC}$ and $M_{SC}$.  Note that the migration rate 

$M_{AC} = m_{AC} N_C$

is the expected number of immigrants in population C that are from population A per generation under the real-world view with time running forward, where $m_{AC}$ is the proportion of immigrants in population C that are from population A.  See Appendix A for definitions of migration rates used in different programs. 

#### Prior on migration rates
The control file variable `migprior = alpha beta` specifies a default gamma prior density with parameters `alpha` and `beta` for all migration rates.  In the example above, both $M_{AC}$ and $M_{SC}$ are assigned the gamma prior $G(2, 200)$, with prior mean $2/200 = 0.01$.
A separate gamma prior can instead be defined for each migration rate (where the migration connection is specified using the source and target populations) and this specification takes precedence.  For example if we specify
```
     migprior = 2 200
    migration = 2
                A C  2 100
                S C
```
the priors are $M_{AC} ~ G(2, 100)$ with the prior mean $0.02$ and $M_{SC} ~ G(2, 200)$ from the default prior specified by `migprior`.

#### Variable migration rates among loci
In this model, the migration rate $M_i$ at each locus *i* varies according to a gamma distribution $M_i \approx G(\alpha_M,\alpha_M/\overline{M})$, with shape parameter $\alpha_M$ while the mean rate $\overline{M}$ is assigned the gamma prior $\overline{M} \approx G(\alpha,\beta)$.  Here $\alpha_M$ is a parameter that characterizes the variation of $M_i$ among loci, with a small $\alpha_M$ (e.g., 0.5 or 1) indicating highly variable rates among loci and a large $\alpha_M$ indicating nearly constant migration rates among loci (when $\alpha_M = \infty$ all loci have the same rate).
```
     migprior = 2 200
    migration = 2
                A C  2 100 5
                S C
```
In the above example, $M_{SC}$ is applied to all loci with the prior $G(2, 200)$, but $M_{AC}$ varies among loci according to the shape parameter 
$\alpha_M = 5$, and the mean rate for all loci $\overline{M}_{AC}$ is assigned the gamma prior $G(2, 100)$.

#### Priors and pseudo-priors
When the gene flow involves ancestral species, migration may become impossible because of changes to the species divergence times in the MCMC.  In the example above, if $\tau_S < \tau_T$ migration from S to C is possible during the time period $(\tau_S, \tau_T)$, but if $\tau_S > \tau_T$ migration from S to C is impossible as the two populations were never contemporary.  Thus the MCMC proposal to change species divergence times $\tau_S$ or $\tau_T$ may cause the migration rate parameter $M_{SC}$ to disappear or reappear.  When $M_{SC}$ is absent from the model (that is, when $\tau_S > \tau_T$), the MCMC algorithm treats it as a pseudo-parameter (written as $M^*_{SC}$) and assigns it a pseudo-prior to facilitate the trans-dimensional move.  The choice of pseudo-priors affects the mixing efficiency of the MCMC, but not the correctness of the algorithm.  For good mixing, one should choose the pseudo-prior to be close to the posterior of the parameter $M_{SC}$ (which can be generated by running short chains).
```
     migprior = 2 200 
    migration = 2
                A C  2 100 
                S C  2 200   100 250
```
In the example above, the prior $M_{SC} \approx G(2, 200)$ is applied when $\tau_S < \tau_T$, and the pseudo-prior $M^*_{SC} \approx G(100, 250)$ with mean 0.4, is applied when $\tau_S > \tau_T$.  Note that in this example the pseudo-prior is far more concentrated than the prior (with shape parameter 100 versus 2).
In order to summarize the posterior from the MCMC samples, only those samples taken when the migration rate is defined should be used (that is, only samples collected when $\tau_S < \tau_T$ in our example).  Samples collected when the migration rate does not exist in the model (that is, when $\tau_S > \tau_T$ in our example) approximate the pseudo-prior.

In the case of four species, the following awk commands may be used to split the samples into two files (this assumes that tau_S and tau_T are in column 10 and 11)
```
awk 'BEGIN {cnt=0} NR>1{if ($10 > $11) {cnt=cnt+1;}}} END {print cnt}' mcmc.txt
# split sample file into ab.txt (samples where t_AB > t_CD) and cd.txt (t_CD > t_AB)
head -n 1 mcmc.txt > cd.txt; awk 'NR>1 {if ($10<$11) print $0}' mcmc.txt >> cd.txt
head -n 1 mcmc.txt > ab.txt; awk 'NR>1 {if ($10>$11) print $0}' mcmc.txt >> ab.txt
```
Use something like this if you understand the syntax. The running mean of the migration rate printed on the monitor during the MCMC run is the posterior mean after the filtering.
Note that this problem of the migration rate appearing and disappearing during the MCMC may exist when a migration event involves ancestral populations, and is not limited to the balanced species tree for four species.  Later we should automate the summary of the MCMC sample.  
Under the model of variable $M$ among loci, all migration rates for all loci may appear and disappear if the migration involves ancestral species.  The following five ways of specifying the model and priors and pseudo-priors are accepted (with $M_{SC}$ as an example)
```
(a)  S C
(b)  S C a_M
(c)  S C a b
(d)  S C a b a_M
(e)  S C a b pseudo_a pseudo_b
(f)  S C a b a_M pseudo_a pseudo_b
```
In options (a) and (b), the default prior specified with migprior is assigned on $M_{SC}$.

## Combined Analyses of Organelle Genomes and Sex Chromosomes
When performing a combined analysis of data from autosomal, mitochondrial, chloroplast or sex chromosomes
one needs to account for the differences of effective population sizes between loci from these
different sources. This can be done using a *heredity multiplier* (referred to as an inheritance scalar by [Hey and Nielsen 2004](https://doi.org/10.1534/genetics.103.024182)). The `heredity` variable in the BPP control file allows for 
such differences. Here are some examples:
```
    heredity = 0       # (0: No variation)
    heredity = 1 4 4   # (1: estimate, & a_gamma b_gamma)
    heredity = 2 heredity.txt      # (2: from file)
```
The specification `heredity = 0` is the default and means that $\theta$ is the same for
all loci. `heredity = 1` or `2` specifies two models that allow $\theta$
to vary among loci, which may be useful for combined analysis of data
from autosomal, mitochondrial, X and Y loci. With such mixed data, the
effective population sizes are different among loci, so that a heredity
multiplier ([Hey and Nielsen 2004](https://doi.org/10.1534/genetics.103.024182)) should be applied. Other
factors such as natural selection may also cause $\theta$ to deviate
from the neutral expectation. BPP implements two options
for this. The first option (`heredity = 1`) is to estimate the
multipliers from the data, using a gamma prior with parameters $\alpha$
and $\beta$ specified by the user. In the example above, a gamma prior
$G(4, 4)$, with mean $4/4 = 1$, is specified for the multiplier for each
locus. The MCMC should then generate a posterior for the multiplier for
each locus. The second option (`heredity = 2`) is for the user to
specify the multipliers in a file, and the multipliers will then be used
as fixed constants in the MCMC run. The file simply contains as many
numerical values as the number of loci, separated by spaces or line
breaks.

| Genome        | Heredity scalar |
|:--------------|:----------------|
| Autosome      | 1               |
| X chromosome  | 0.75            |
| Y chromosome  | 0.25            |
| Mitochondrial | 0.25            |

Table. Examples of Common Heredity Scalars

Note.--- The effect of the locus-specific mutation rates and the
locus-specific heredity multipliers are different. A locus rate is used
to multiply all $\theta$s and $\tau$s for the locus, while a heredity
multiplier is used to multiply all $\theta$ parameters for the locus but
not the $\tau$s. Nevertheless, those parameters are quite likely to be
strongly correlated, especially when the species tree is small.

## Methods of Analysis

Here we describe the specifics of the four different types of analyses
that BPP can perform and illustrate them using several of the example
control files provided in the `examples` subdirectory of the BPP
distribution. Briefly, BPP is capable of analyzing sequence data using 4
distinct phylogeographical models. The 4 models are determined by the
combinations of two possible settings for each of the two Boolean
control file variables `speciesdelimitation` (0 = no species
delimitation, 1 = species delimitation) and `speciestree` (0 = fixed
species tree, 1 = inferred species tree). Therefore, analysis A00 infers
the biogeographical parameters $\theta$ and $\tau$ on a fixed (user
specified) tree with a fixed (delimited) number of species/populations
as described in ([Rannala and Yang 2003](https://doi.org/10.1093/genetics/164.4.1645). Analysis A10 delimits the number of
species and infers biogeographical parameters using a fixed species tree
(guide tree) as described in [Yang and Rannala 2010](https://doi.org/10.1073/pnas.0913022107). Analysis A01 infers the
species tree and biogeographical parameters using a fixed number of
species (delimitation) as described in [Rannala and Yang 2017](https://doi.org/10.1093/sysbio/syw119). Analysis A11
jointly infers the species tree, species delimitation and
biogeographical parameters as described in [Yang and Rannala 2014](https://doi.org/10.1093/molbev/msu279). 
Under Model A00 it is also possible to include a model of either instantaneous
introgression or ongoing gene flow -- examples of such analyses are
given in [Introgression and Migration Models](#introgression-and-migration-models). For each of these 4 major types of analyses
the control file specifications will differ as well as the form of
output printed to the screen and output files. Therefore, we will
include two subsections in our discussion of each analysis method: input
and output. Here we only consider the control file variables specific to
the methods of analysis under consideration. For a general discussion of
other control file variables see the example control file discussed in [Control File](#control-file).
Note that for most analyses the output printed to screen during the run
and that printed to the specified output file will be identical.

It is important to note that only analysis A00 applies within-model
inference. There is one specified model in that case (the MSC with a
fixed species tree) and the parameters are well defined, with the
objective being to estimate those parameters. In contrast, analyses A01,
A10, and A11 all apply trans-model inference. They move between
different models, and the main objective is to calculate the posterior
probabilities of the models. Each is an instance of the MSC model, but
the species delimitation (the number and nature of the species) and/or
the species phylogeny may differ between models. In analyses A01, A10,
and A11, the prior specified using control file variable `thetaprior`
applies to all $\theta$ parameters in all models. Similarly, all MSC
models specifying two or more delimited species have a parameter
$\tau_0$ (the divergence time of the root), and these parameters (for
different models) are assigned the same prior, specified by the
`tauprior` control file variable.

### A00: Demography and Divergence Time Estimation 

Analysis A00`(speciesdelimitation = 0, speciestree = 0)`, requires a
user-specified species tree topology, given by variable `species&tree`
in the control file, and infers the $\theta$ and $\tau$ parameters on
this fixed tree. The theory underlying estimation of species divergence
times and population sizes ($\tau$s and $\theta$s) under the MSC model
when the species phylogeny is given is described in ([Rannala and Yang 2003](https://doi.org/10.1093/genetics/164.4.1645). Note
that if there are $s$ species in the species tree, the model will
minimally contain the following parameters: $s – 1$ species divergence
times ($\tau$s) and $s – 1$ ancestral species $\theta$s. A $\theta$
parameter will also exist for each contemporary species with at least
two sequences present at some locus. If the data are for a single
species/population, the model will contain one parameter only, $\theta$
for that species. If there is only one species, the MSC model becomes
the single-population Kingman coalescent. We will provide
examples of both a single population analysis (using the control file
`yu2001.bpp.ctl`) and a multi-species analysis (using the control file
`A00.bpp.ctl`).

**How many $\theta$ and $\tau$ parameters exist?**  
The number of parameters that BPP includes in the MSC or
MSC-I/MSC-M models depends on the data configuration. The program always
includes $\theta$ and $\tau$ parameters for each internal node
(ancestral species) on the species tree, but includes a $\theta$ for an
extant species (tip) if and only if that species has at least 2
sequences at some loci. If two or more sequences for a species are
specified in the control file but there are at most 0 or 1 sequence per
locus for that species in the sequence data file, $\theta$ for that
species will not be identifiable or estimable. In that case, the
posterior for the parameter will be the prior. Nevertheless, other
results (such as the posterior distribution for other parameters or
posterior probabilities for species trees and species delimitations)
will still be correct. The same applies to other more complex cases of
missing data and parameter non-identifiability. As an example, suppose
the species tree and the numbers of sequences (for species A, B, C and D from left to right) for two kinds of loci are
as follows:

((A,B), (C,D))  
locus configuration 1: 1 1 0 0  
locus configuration 2: 2 0 0 1

In this case, $\theta_A$, $\theta_{AB}$, $\theta_{ABCD}$, $\tau_{AB}$,
and $\tau_{ABCD}$ are identifiable while $\theta_B$, $\theta_C$,
$\theta_D$, $\theta_{CD}$ and $\tau_{CD}$ are not. It is impossible to
estimate $\theta_D$ , $\theta_{CD}$ and $\tau_{CD}$ as no data are
available from species C.

#### Input A00 (single population)

Our first example is a dataset of human sequence data generated by
[Yu et al 2001](https://doi.org/10.1093/oxfordjournals.molbev.a003795). 
The control file is found in the BPP distribution at location:
```
    examples/yu2001/yu2001.bpp.ctl
```
The data comprise 61 phased sequences for a single locus that will be
analyzed to estimate the single parameter $\theta$. There is no need to
specify an Imap file, or tag the sequence names in the sequence file
(`yu2001.txt`); the sequence names are read and then ignored. Multiple
loci may be included in the sequence file. The contents of the control
file `yu2001.bpp.ctl` are shown below:
```
    seed = -1

    seqfile = yu2001.txt
    outfile = out.txt
    mcmcfile = mcmc.txt

    # fixed species delimitation and species tree
    speciesdelimitation = 0 
    speciestree = 0

    species&tree = 1  H
                      61  

    # 0: no data (prior); 1:seq like
    usedata = 1

    # number of data sets in seqfile
    nloci = 1    

    # remove sites with ambiguity data (1:yes, 0:no)?
    cleandata = 0    

    # gamma(a, b) for theta
    thetaprior = gamma 2 2000

    # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
    finetune = 1: 2 0.00001 0.0001  0.0005 0.5 0.2 1.0  

    # MCMC samples, locusrate, heredityscalars, genetrees, locussubstitution
    print = 1 0 0 0  

    burnin = 4000
    sampfreq = 2
    nsample = 10000
```
The two lines:
```
    speciesdelimitation = 0 
    speciestree = 0
```
specify that the species delimitation (number of species) is fixed as
well as the species tree topology. In this case, there is no species
tree and only a single species exists. There is no need for a species
tree, so the block for specifying species names and species tree looks
like this:
```
    species&tree = 1  H
                      61
```
The line:
```
    # gamma(a, b) for theta
    thetaprior = gamma 2 2000
```
specifies a gamma distribution for the prior on $\theta$ with parameters
$\alpha = 2$ and $\beta = 2000$ which has mean $0.001$ and variance
$5 \times 10^{-7}$.

#### Output A00 (single population)

When BPP is run using the above control file a summary of the progress
of the MCMC analysis is printed to screen. Here, we briefly explain what
this output means. The different analyses presented in this chapter all
produce a similar form of output when the program is running, but the
number of columns and their specific content varies somewhat. Omitting
some mostly irrelevant output related to the adjustment of proposal
moves, the output to screen when this control file is run appears as
follows:
```{ .yaml .no-copy }
         |  Acceptance proportions  |
    Prgs | Gage Gspr thet  tau  mix | mthet1       log-PG         log-L
    -------------------------------------------------------------------
     -15%  0.71 0.18 0.30 0.00 0.38   0.0004    467.11472  -12721.22305  0:01
     (some output omitted here)
       5%  0.71 0.29 0.33 0.00 0.24   0.0004    428.40742  -12721.09413  0:05
      10%  0.71 0.29 0.31 0.00 0.25   0.0003    480.92848  -12720.75185  0:06
      15%  0.71 0.29 0.30 0.00 0.24   0.0003    476.37432  -12720.65199  0:08
      20%  0.71 0.29 0.31 0.00 0.23   0.0003    441.36960  -12720.92326  0:09
      25%  0.71 0.29 0.32 0.00 0.23   0.0003    443.51740  -12720.85802  0:10
      30%  0.71 0.29 0.31 0.00 0.23   0.0003    441.93616  -12720.86033  0:11
      35%  0.71 0.29 0.31 0.00 0.23   0.0003    495.50716  -12720.84565  0:12
      40%  0.71 0.29 0.31 0.00 0.23   0.0003    452.82049  -12720.87810  0:13
      45%  0.71 0.29 0.31 0.00 0.23   0.0003    461.93346  -12720.94582  0:15
      50%  0.71 0.29 0.32 0.00 0.23   0.0004    463.85955  -12721.01005  0:16
      55%  0.71 0.29 0.32 0.00 0.23   0.0004    463.42744  -12720.98006  0:17
      60%  0.71 0.29 0.32 0.00 0.23   0.0004    478.52797  -12721.03353  0:18
      65%  0.71 0.29 0.32 0.00 0.23   0.0004    499.97705  -12720.96168  0:19
      70%  0.71 0.29 0.32 0.00 0.23   0.0004    455.32489  -12720.93611  0:20
      75%  0.71 0.29 0.32 0.00 0.23   0.0004    439.21544  -12720.99190  0:22
      80%  0.71 0.29 0.32 0.00 0.23   0.0004    461.56954  -12721.00180  0:23
      85%  0.71 0.29 0.32 0.00 0.23   0.0004    469.53443  -12720.99823  0:24
      90%  0.71 0.29 0.32 0.00 0.23   0.0004    455.60303  -12720.99213  0:25
      95%  0.71 0.29 0.32 0.00 0.23   0.0004    454.07141  -12721.00745  0:26
     100%  0.71 0.29 0.32 0.00 0.23   0.0004    461.01044  -12720.97549  0:27

    0:27 spent in MCMC

              theta_1H  lnL
    mean      0.000352  -12720.995649
    median    0.000337  -12720.726000
    S.D       0.000116  2.936769
    min       0.000089  -12735.535000
    max       0.001293  -12712.774000
    2.5%      0.000170  -12727.454000
    97.5%     0.000617  -12716.086000
    2.5%HPD   0.000145  -12726.920000
    97.5%HPD  0.000574  -12715.718000
    ESS*      757.838791  1180.434956
    Eff*      0.075784  0.118043
```
Note that in these examples we are using a random seed from the computer clock so your results will differ slightly from those shown. The line:
```{ .yaml .no-copy }
         |  Acceptance proportions  |
    Prgs | Gage Gspr thet  tau  mix | mthet1       log-PG         log-L
    -------------------------------------------------------------------
```
is a header that explains the content of each column printed during the
run. The `Prgs` (progress) column will appear in all analyses and
indicates the percentage of the MCMC iterations that have been
completed. If this number is negative it indicates that the program is
still running the burn-in iterations. For example, `-15%` means that 15
percent of the burn-in iterations remain to be completed. The next 5
columns are the current acceptance proportions for different parameter
proposals in the MCMC. The proposals are defined as follows:

-   `Gage`: proposal to change ages of nodes on gene trees

-   `Gspr`: proposal to change gene tree topology using SPR move

-   `thet`: proposal to change $\theta$ parameter

-   `tau`: proposal to change $\tau$ parameter

-   `mix`: mixing step proposal

In the output the acceptance proportions are stable. The fifth column
(the $\tau$ proposal acceptance proportion) is zero; this is expected
because there are no species divergence times ($\tau$s) in the model
when only a single species exists. Column 7, labeled `mthet1` gives the
current average (mean) value of theta in the MCMC. Columns 8 and 9 are
the log probability of the coalescent model (`log-PG`) and the
log-likelihood of the sequence data on the gene tree (`log-L`). Check
that the acceptance proportions are in a reasonable range (about 20% to
40% for continuous parameters such as $\theta$). If you need to
terminate the program before it finishes, you can do so using the key
combination `Ctrl-C` (e.g., hold down the Control key and then press the
C key). The final screen output, printed after the MCMC is completed,
summarizes the posterior parameter estimates for $\theta$ (`theta_1H`)
in the second column. The posterior mean in the first row, for example,
is $0.000352$ and the lower and upper limits of the highest posterior
density (HPD) set of values are in rows 8 and 9, respectively. In this
case, the HPD credibility interval for $\theta$ is
$(0.000145,0.000574)$.

#### Input A00 (multiple populations)

Our second example is a dataset of frog sequence data generated by
[Zhou et al 2012](https://doi.org/10.1111/j.1365-294X.2011.05411.x). 
This same dataset will also be used to illustrate analyses
A01, A10 and A11. The control file is found in the BPP distribution at
location:
```
    examples/frogs/A00.bpp.ctl
```
The data comprise 5 loci for 4 populations, with the number of sequences
at each locus varying between 21 and 30. We will be estimating $\theta$s
for each of the 4 contemporary and 3 ancestral populations as well as
the 3 species divergence times ($\tau$s). The contents of the control
file `A00.bpp.ctl` are shown below:
```
    seed =  -1

    seqfile = frogs.txt
    Imapfile = frogs.Imap.txt
    outfile = out.txt
    mcmcfile = mcmc.txt

    # fixed number of species/populations
    speciesdelimitation = 0 

    # fixed species tree
    speciestree = 0

    species&tree = 4  K  C  L  H
                        9  7 14  2
                     (((K, C), L), H);

    # unphased data for all 4 populations
    phase =   1  1  1  1

    # use sequence likelihood                  
    usedata = 1

    nloci = 5  

    # do not remove sites with ambiguity data
    cleandata = 0    

    # gamma(a, b) for theta (estimate theta)
    thetaprior = gamma 2 2000 

    # gamma(a, b) for root tau & Dirichlet(a) for other tau's
    tauprior = gamma 2 1000 

    # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
    finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  

    # MCMCsamples, locusrate, heredityscalars, genetrees, substitutionparams
    print = 1 0 0 0

    burnin = 8000
    sampfreq = 2
    nsample = 100000
```
As in the previous analysis, the two lines:
```
    speciesdelimitation = 0 
    speciestree = 0
```
specify that the species delimitation (number of species) is fixed as
well as the species tree topology. In this case, a species tree is
specified along with population labels. The block for specifying species
names and the species tree looks like this:
```
    species&tree = 4  K  C  L  H
                        9  7 14  2
                     (((K, C), L), H);
```
The first line above lists the number of species (populations) which is
4 in this example. This is followed by the label for each species (these
should match the labels in the map file):
```
    K  C  L  H
```
The next line specifies the maximum number of sequences per species in
the same order as the species labels:
```
    9  7 14  2
```
the actual number of sequences at any locus must be less than this
value. The third line is a rooted phylogenetic tree topology specified
in Newick format:
```
    (((K, C), L), H);
```
The next line specifies that the loci for each of the 4 populations are
unphased:
```
    phase = 1 1 1 1
```
The specified phase variables should include a 0 or 1 (Boolean) entry
for each population, where the populations are assumed to be in the same
order as specified above (in the `species&tree` block) and the Boolean
variable indicates that the sequence data for a population are either
unphased (1) or phased (0). Ambiguity characters are used to represent
genotypes for unphased sequences (see [BPP control file variables](#bpp-control-file-variables)). 
In the example data the sequences
are unphased for all 4 populations. Currently, it is not possible to
have partially phased data (e.g., only some loci unphased in a
particular population). The program will infer the phase as part of the
analysis.

Unlike the earlier single population analysis, in which we estimated a
single $\theta$ parameter, in the multipopulation analysis we will be
estimating multiple $\theta$s (one for each population) and divergence
times ($\tau$s) between populations. Thus, we have prior distributions
for both $\theta$ and $\tau_0$. The lines below specify a gamma
distribution as the prior for $\theta$ with parameters $\alpha = 2$ and
$\beta = 2000$ (the prior mean and variance of $\theta$ are $0.001$ and
$0.0000005$):
```
    # gamma(a, b) for theta
    thetaprior = gamma 2 2000 
```
The lines:
```
    # gamma(a, b) for root tau & Dirichlet(a) for other tau's
    tauprior = gamma 1 1000 
```
specify a gamma distribution as a prior for the root age and
the remaining $\tau$s have a uniform Dirichlet distribution conditional
on the root age (the prior mean and variance of the root age are $0.001$
and $0.000001$,respectively):

#### Output A00 (multiple populations)

When BPP is run using the above control file a summary of the progress
of the MCMC is again printed to screen as follows:
```{ .yaml .no-copy }
         |  Acceptance proportions  |
    Prgs | Gage Gspr thet  tau  mix | mthet1 mthet2 mthet3   mtau1  mtau2  mtau3       log-PG        log-L
    ------------------------------------------------------------------------------------------------------
      -3%  0.65 0.10 0.43 0.04 0.12   0.0031 0.0095 0.0068  0.0021 0.0015 0.0014   1204.63468  -4458.38814
     (some output omitted here)
       5%  0.64 0.26 0.30 0.28 0.27   0.0035 0.0096 0.0067  0.0019 0.0018 0.0017   1225.51650  -4443.85539  1:09
      10%  0.65 0.27 0.30 0.28 0.28   0.0035 0.0095 0.0068  0.0019 0.0018 0.0018   1202.62497  -4445.50020  1:47
      15%  0.64 0.27 0.30 0.28 0.28   0.0034 0.0095 0.0068  0.0019 0.0018 0.0017   1211.95125  -4446.35694  2:25
      20%  0.64 0.26 0.31 0.28 0.28   0.0034 0.0095 0.0067  0.0018 0.0018 0.0017   1211.09206  -4445.74292  3:05
      25%  0.64 0.26 0.31 0.28 0.28   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1222.97603  -4446.12020  3:43
      30%  0.64 0.26 0.31 0.28 0.28   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1216.60217  -4445.55424  4:21
      35%  0.64 0.26 0.31 0.28 0.28   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1260.59722  -4445.22113  5:00
      40%  0.64 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1269.45020  -4445.20340  5:38
      45%  0.65 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1243.07668  -4445.00207  6:16
      50%  0.65 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1229.59858  -4445.04036  6:55
      55%  0.64 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1209.31273  -4445.43161  7:33
      60%  0.64 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1206.22162  -4445.53994  8:11
      65%  0.64 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1190.78837  -4445.21618  8:50
      70%  0.64 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1203.82684  -4445.04124  9:30
      75%  0.64 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1186.32176  -4445.00578  10:09
      80%  0.64 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1222.75718  -4444.94315  10:47
      85%  0.64 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1228.39106  -4444.68671  11:26
      90%  0.64 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1228.72385  -4444.54844  12:04
      95%  0.64 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1223.25921  -4444.25455  12:43
     100%  0.64 0.26 0.31 0.28 0.27   0.0034 0.0096 0.0067  0.0018 0.0018 0.0017   1199.52731  -4444.17857  13:21
    }
    13:21 spent in MCMC

              theta_1K  theta_2C  theta_3L  theta_4H  theta_5KCLH theta_6KCL theta_7KC tau_5KCLH tau_6KCL tau_7KC lnL
    mean      0.003422  0.009588  0.006683  0.003031  0.003720  0.001577  0.001648  0.001823  0.001787  0.001710  -4444.16
    median    0.003351  0.009492  0.006599  0.002932  0.003642  0.001412  0.001463  0.001809  0.001773  0.001700  -4443.69
    S.D       0.000747  0.001461  0.001118  0.000819  0.000944  0.000899  0.000954  0.000234  0.000230  0.000238  11.615936
    min       0.001309  0.004913  0.003258  0.000866  0.001131  0.000029  0.000039  0.001047  0.001043  0.000772  -4496.72
    max       0.008147  0.017569  0.012778  0.009724  0.010355  0.008461  0.008987  0.003140  0.003113  0.002847  -4403.99
    2.5%      0.002171  0.006999  0.004743  0.001729  0.002100  0.000325  0.000340  0.001404  0.001375  0.001271  -4468.01
    97.5%     0.005084  0.012720  0.009104  0.004919  0.005780  0.003751  0.003991  0.002319  0.002280  0.002209  -4422.74
    2.5%HPD   0.002071  0.006802  0.004595  0.001566  0.001949  0.000142  0.000164  0.001396  0.001342  0.001252  -4467.06
    97.5%HPD  0.004928  0.012469  0.008903  0.004650  0.005573  0.003312  0.003512  0.002308  0.002237  0.002185  -4421.90
    ESS*      2721.25   8954.56   5279.29   11762.38  1199.63   4543.27   3299.32   850.33    970.61    932.00    420.40
    Eff*      0.027213  0.089546  0.052793  0.117624  0.011996  0.045433  0.032993  0.008503  0.009706  0.009320  0.004204
```
There are again 5 acceptance proportions for the 5 parameter proposals,
but now the current mean values of 6 parameters are printed. The current
mean values of $\theta$ for nodes 1, 2 and 3 (`mthet1 mthet2 mthet3`)
are printed -- which are $\theta$s of the contemporary species K, C and
L, and those of the 3 $\tau$s (`mtau1 mtau2 mtau3`) which are the
divergence times for ancestral species KCLH, KCL and KC, respectively.
If there are many populations/species only the first few $\tau$s will be
printed to screen but all the $\tau$s will be available in the summary
created at the end of the run (see below). When the run finishes, a
final block is written summarizing the results. The first two rows give
the mean and median of the posterior distribution for each parameter.
This is followed by summaries of the statistical uncertainty: the
standard deviation (S.D.) of the posterior distribution; the minimum
(min) and maximum (max) values; the lower and upper bounds of the 95%
[Credible Interval](https://en.wikipedia.org/wiki/Credible_interval)
(2.5% and 97.5%); the lower and upper bounds of the [Highest Posterior
Density (HPD) Interval](https://en.wikipedia.org/wiki/Credible_interval)
(2.5%HPD and 97.5%HPD); the Effective Sample Size (ESS\*); and the
Efficiency (Eff\*). Larger values for ESS\* and Eff\* indicate better
mixing and more reliable estimates. The above information is also
printed (along with other technical details of the run) to the output
file specified by the `outfile` variable in the control file (in our
example, `out.txt`). A tree file in NEXUS/Newick format is also created
and placed in the file named FigTree.tre which is formatted for
viewing/printing using the
[Figtree](http://tree.bio.ed.ac.uk/software/figtree/) program

#### MCMC output file

A file will be produced (with a name specified by the control file
variable `mcmcfile`) that contains the mcmc samples from the run for all
the $\theta$ and $\tau$ parameters. When running under analysis method
A00 the contents of this file can be analyzed using a program such as
[Tracer]() to visually examine the convergence and mixing of the MCMC
run. However, the other three analysis methods (A01, A10 and A11) are
trans-model analyses in which the number of parameters and or/the
meaning of the parameters changes as the chain runs so it is not correct
to examine the trace plot without conditioning on a particular model.
An example of 5 lines from the `mcmc.txt` file produced by the above
example are given below:
```{ .yaml .no-copy }
Gen	theta_1K	theta_2C	theta_3L	theta_4H	theta_5KCLH	theta_6KCL	theta_7KC	tau_5KCLH	tau_6KCL	tau_7KC	lnL
2	0.002848	0.010633	0.006000	0.002231	0.003697	0.001594	0.000704	0.001754	0.001717	0.001708	-4466.479
4	0.002848	0.011774	0.004789	0.002231	0.003697	0.002134	0.000798	0.001754	0.001717	0.001708	-4457.097
6	0.002848	0.010660	0.005826	0.002231	0.003697	0.001250	0.000799	0.001754	0.001717	0.001705	-4448.002
8	0.002848	0.012103	0.005826	0.002231	0.002615	0.001250	0.000803	0.001754	0.001717	0.001696	-4455.900
```

### A01: Species Tree Estimation

Analysis A01 `(speciedelimitation = 0, speciestree = 1)` assumes that the species delimitation is fixed
(e.g., species assignments are as specified in the Imap file) but that the species tree is unknown. The
program estimated the species tree topology in addition to the $\theta$s and $\tau$s.

#### Input A01
Our example uses the dataset of frog sequence data generated by [Zhou et al 2012](https://doi.org/10.1111/j.1365-294X.2011.05411.x)
that we considered previously. The contents of the control file `A01.bpp.ctl` are shown below:
```
  seed =  -1
  seqfile = frogs.txt
  Imapfile = frogs.Imap.txt
  outfile = out.txt
  mcmcfile = mcmc.txt

  speciesdelimitation = 0 * fixed species tree
  speciestree = 1  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

  speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 4  K  C  L  H
                    9  7 14  2
                 ((K, C), (L, H));
        phase =   1  1  1  1
   
  usedata = 1  * 0: no data (prior); 1:seq like
  nloci = 5  * number of data sets in seqfile

  cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

  thetaprior = gamma 2 2000 # gamma(a, b) for theta (estimate theta)
  tauprior = gamma 2 1000 # gamma(a, b) for root tau & Dirichlet(a) for other tau's

  finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

  print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
  burnin = 8000
  sampfreq = 2
  nsample = 100000
```
The two lines:
```
  speciesdelimitation = 0 * fixed species delimitation
  speciestree = 1  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio
```
now specify that the species delimitation is fixed and the species tree topology is being estimated
`speciestree = 1`. The 3 parameters following the `speciestree` variable specify the probability that the
*snakes and ladders* (SNL)  move is used, rather than the subtree-pruning-regrafting (SPR) move to modify the species tree
and the expand/shrink ratios for the SNL move; these proportions can be adjusted on the interval $(0,1)$
and will affect mixing of the MCMC. For more details, see the description of the control file variable [14 speciestree](#14-speciestree).
If no parameters are specified after `speciestree = 1` then the SPR move is used exclusively. See 
[Rannala and Yang 2017](https://doi.org/10.1093/sysbio/syw119) for a description of the proposals available
to change species trees in the MCMC.
The line:
```
 speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted
```
specifies the prior distribution on the species tree topologies. The specification `speciesmodelprior = 1` uses
a prior that gives equal probability to all rooted trees apriori. See the description of the control file
variable [15 speciesmodelprior](#15-speciesmodelprior) for more details. The available priors are described in 
[Yang and Rannala 2014](https://doi.org/10.1093/molbev/msu279).

#### Output A01
When BPP is run using the above control file a summary of the progress of the MCMC is again printed to screen as follows:
```{ .yaml .no-copy }

     |         Acceptance proportions         |
Prgs | Gage Gspr thet  tau  mix   Sspr   Ssnl | mthet1   mtau1       log-PG        log-L
---------------------------------------------------------------------------------------
(some output omitted here)
   5%  0.64 0.26 0.31 0.28 0.27 0.0980 0.0000   0.0037  0.0019   1213.12189  -4448.22572  1:00
  10%  0.64 0.26 0.30 0.28 0.27 0.0865 0.0000   0.0037  0.0019   1212.13271  -4447.30801  1:37
  15%  0.64 0.26 0.30 0.27 0.27 0.0885 0.0000   0.0037  0.0019   1205.76722  -4446.44035  2:20
  20%  0.64 0.26 0.30 0.27 0.27 0.0867 0.0000   0.0037  0.0019   1209.18617  -4445.20436  3:03
  25%  0.64 0.26 0.31 0.27 0.27 0.0919 0.0000   0.0038  0.0019   1189.24742  -4444.02665  3:45
  30%  0.65 0.26 0.31 0.27 0.27 0.0959 0.0000   0.0038  0.0018   1218.30083  -4443.50919  4:28
  35%  0.65 0.26 0.31 0.27 0.27 0.0908 0.0000   0.0038  0.0019   1271.24840  -4443.30946  5:13
  40%  0.65 0.26 0.31 0.27 0.27 0.0903 0.0000   0.0038  0.0018   1221.51424  -4442.93455  5:56
  45%  0.65 0.26 0.31 0.28 0.27 0.0865 0.0000   0.0037  0.0018   1179.64791  -4442.58322  6:39
  50%  0.65 0.26 0.31 0.28 0.27 0.0883 0.0000   0.0038  0.0018   1220.10587  -4442.33168  7:21
  55%  0.65 0.26 0.31 0.28 0.27 0.0882 0.0000   0.0037  0.0018   1230.47093  -4442.53042  8:03
  60%  0.65 0.26 0.31 0.28 0.27 0.0900 0.0000   0.0037  0.0019   1210.69277  -4443.13901  8:46
  65%  0.65 0.26 0.31 0.28 0.27 0.0906 0.0000   0.0037  0.0019   1259.02426  -4443.42635  9:29
  70%  0.65 0.26 0.31 0.28 0.27 0.0907 0.0000   0.0037  0.0019   1227.92007  -4443.74417  10:11
  75%  0.65 0.26 0.31 0.28 0.27 0.0907 0.0000   0.0037  0.0019   1254.82843  -4443.96093  10:55
  80%  0.65 0.26 0.31 0.28 0.27 0.0912 0.0000   0.0037  0.0019   1211.86344  -4443.93739  11:40
  85%  0.65 0.26 0.31 0.28 0.27 0.0922 0.0000   0.0037  0.0019   1205.06522  -4443.99758  12:23
  90%  0.65 0.26 0.31 0.28 0.27 0.0924 0.0000   0.0037  0.0019   1191.95318  -4444.06010  13:05
  95%  0.65 0.26 0.31 0.28 0.27 0.0918 0.0000   0.0037  0.0019   1213.34484  -4444.05286  13:47
 100%  0.65 0.26 0.31 0.28 0.27 0.0916 0.0000   0.0037  0.0019   1224.45872  -4443.84148  14:30

14:30 spent in MCMC

Species in order:
   1. K
   2. C
   3. L
   4. H

(A) Best trees in the sample (15 distinct trees in all)
    33978  0.33978  0.33978 ((C, (H, L)), K);
    17928  0.17928  0.51905 (C, ((H, L), K));
     8961  0.08961  0.60866 ((C, K), (H, L));
     8236  0.08236  0.69102 (C, ((H, K), L));
     6896  0.06896  0.75998 (((C, H), L), K);
     4548  0.04548  0.80546 (((C, L), H), K);
     4066  0.04066  0.84612 (C, (H, (K, L)));
     3037  0.03037  0.87649 ((C, (H, K)), L);
     2775  0.02775  0.90424 (((C, K), H), L);
     2245  0.02245  0.92669 (((C, H), K), L);
     2131  0.02131  0.94800 ((C, L), (H, K));
     1786  0.01786  0.96586 (((C, K), L), H);
     1423  0.01423  0.98009 (((C, L), K), H);
     1094  0.01094  0.99103 ((C, H), (K, L));
      897  0.00897  1.00000 ((C, (K, L)), H);

(B) Best splits in the sample of trees (10 splits in all)
 60867 0.608664  0011
 45422 0.454215  0111
 30230 0.302297  1011
 13522 0.135219  1100
 13404 0.134039  1001
 10235 0.102349  0101
  8102 0.081019  0110
  8057 0.080569  1101
  6057 0.060569  1010
  4106 0.041060  1110

(C) Majority-rule consensus tree
(K, C, (L, H) #0.608664);

(D) Best tree (or trees from the mastertree file) with support values
((C, (H, L) #0.608664) #0.454215, K);   [P = 0.339777]
```
The first part of this output are the current acceptance proportions and running means
for different parameters as the MCMC was run. 
Columns 2 to 6 are the five acceptance proportions for the conventional MCMC moves, as discussed above. 
Columns 7 and 8 are the acceptance proportions for moves that
change the species tree topology, both are much lower than the
other parameter proposals, this is typical and it is often not possible
to improve these acceptance proportions. 
Columns 9 and 10 are posterior means of $\theta$ and $\tau$ for the root population (-1 is printed in column 9 if $\theta$s are integrated out). The last two numbers are the log MSC gene-tree density ([Rannala and Yang 2003](https://doi.org/10.1093/genetics/164.4.1645)) and the average log sequence likelihood (Felsenstein 1981).

Next, there are 4 sections: A, B, C and D summarizing the results.
Section (A) lists the species trees in decreasing order of posterior
probabilities. From this you can easily identify the 95% or 99%
credibility set of species trees. Section (B) lists the splits (or
bipartitions) and their posterior probabilities. The splits here take
into account the location of the root, and may be different from the
splits for unrooted trees. Section (C) prints the majority-rule
consensus tree, with posterior probabilities for nodes.

#### MCMC output file
The MCMC sample of species trees is collected in the file mcmc.txt.
Below are five lines from that file. The numbers after `:` are the branch
lengths ($\tau$s), while those after `#` are the $\theta$s. 
```{ .yaml .no-copy }
(K #0.002983: 0.001998, ((C #0.009671: 0.001684, H #0.003031: 0.001684) #0.001835: 0.000250, L #0.006266: 0.001934) #0.000560: 0.000065) #0.002770;
(K #0.003247: 0.002176, ((C #0.010528: 0.001871, H #0.004726: 0.001871) #0.002902: 0.000147, L #0.006615: 0.002017) #0.000950: 0.000158) #0.003016;
(K #0.003234: 0.002167, ((C #0.010486: 0.001843, H #0.003294: 0.001843) #0.002186: 0.000147, L #0.007822: 0.001989) #0.001087: 0.000178) #0.003003;
(K #0.002677: 0.001794, ((C #0.010973: 0.001516, H #0.001906: 0.001516) #0.001821: 0.000093, L #0.005827: 0.001609) #0.001261: 0.000185) #0.002486;
(K #0.002677: 0.001794, ((C #0.009869: 0.001516, H #0.003071: 0.001516) #0.001821: 0.000109, L #0.005788: 0.001625) #0.002495: 0.000169) #0.002486;
```

### A10: Species Delimitation
Analysis A10 `(speciedelimitation = 1, speciestree = 0)` assumes that the species delimitation is unknown
(e.g., population assignments as specified in the Imap file may not represent distinct species) but that the species *guide tree* is fixed (known). The
program calculates the posterior probabilities of different possible species delimitations as well as of $\theta$s and $\tau$s.

#### Input A10
Our example uses the dataset of frog sequence data generated by [Zhou et al 2012](https://doi.org/10.1111/j.1365-294X.2011.05411.x)
that we considered previously. The contents of the control file `A10.bpp.ctl` are shown below:
```
	seed =  -1
	seqfile = frogs.txt
	Imapfile = frogs.Imap.txt
	outfile = out.txt
	mcmcfile = mcmc.txt

	speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
	speciestree = 0        * species tree NNI/SPR

	speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

	species&tree = 4  K  C  L  H
                      9  7 14  2
                 ((K, C), (L, H));
	phase =   1  1  1  1
                  
	usedata = 1  * 0: no data (prior); 1:seq like
	nloci = 1  * number of data sets in seqfile

	cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

	thetaprior = gamma 2 2000 # gamma(a, b) for theta (estimate theta)
	tauprior = gamma 2 1000 # gamma(a, b) for root tau & Dirichlet(a) for other tau's

	finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

	print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
	burnin = 8000
	sampfreq = 2
	nsample = 100000
```
the lines:
```
	speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
	speciestree = 0        * species tree NNI/SPR
```
specify that the species tree is fixed `speciestree = 0` (a guide tree is used as specified in the control file)
and the rjMCMC algorithm 1, with $\alpha = 2$ and $m = 1$ in equations 6 and 7 of 
[Yang and Rannala 2010](https://doi.org/10.1073/pnas.0913022107) is used for species
delimitation. Another possible specification for delimitation would be:
```
speciesdelimitation = 1 0 2
```
which would specify rjMCMC algorithm 0 with $\epsilon = 2$ in
equations 3 and 4 of [Yang and Rannala 2010](https://doi.org/10.1073/pnas.0913022107).
The two algorithms in theory should produce identical results. 
The line:
```
speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted
```
specifies the prior distribution on topologies (see above) Prior 0 means
equal probabilities for labeled histories (which are rooted trees with
internal nodes ordered by their age). This is the prior used by
[Yang and Rannala 2010](https://doi.org/10.1073/pnas.0913022107).
Prior 1 means equal probabilities for rooted trees (now the default). 
The prior with user specified probabilities
for nodes described by [Rannala and Yang 2013](https://doi.org/10.1534/genetics.112.149039) was removed from BPP after version 2.3.
We also reduced the number of loci to 1 (`nloci = 1`) in this control file because otherwise species
delimitation model 5 (three species) has a posterior probability of 1.

#### Output A10
The species delimitation models that can be generated from the fixed guide tree will be
listed in the output, together with their prior probabilities calculated by
BPP. (As a check, if you use `usedata = 0`, the MCMC
should be sampling from this prior distribution.) The species
delimitation model (for 4 populations as in this data set) is represented using three 0-1 flags for the three
interior nodes 5, 6, 7 in the guide tree, with 0 for 'collapsed' and
1 for 'resolved'. Note that the tips in the guide tree are numbered 1,
2, $\cdots, s$ for $s$ potential species, while the interior (ancestral)
nodes are numbered $s + 1, s + 2, \cdots, 2s – 1$, with $s + 1$ to be
the root of the guide tree. The numbering is through a tree-traversal
algorithm, fixed by the program. This same order is used to specify the
divergence time parameters ($\tau$s), so you can work out the order by
looking at the list of nodes in the screen output (look at the
"population by population table", "# species divergence times in the
order:", etc.).

When BPP is run using the above control file a summary of the progress of the MCMC is again printed to screen as follows:
```{ .yaml .no-copy }
     |     Acceptance proportions      |
Prgs | Gage Gspr thet  tau  mix     rj | np del        mldp mthet1   mtau1       log-PG        log-L
----------------------------------------------------------------------------------------------------
(some output omitted here)
   5%  0.62 0.25 0.30 0.33 0.29 0.0061   10 111 P[5]=0.9726 0.0026  0.0009   217.01433   -972.53408  0:08
  10%  0.63 0.25 0.30 0.34 0.29 0.0037   10 111 P[5]=0.9845 0.0026  0.0009   229.43024   -972.60186  0:12
  15%  0.63 0.25 0.30 0.34 0.29 0.0033   10 111 P[5]=0.9886 0.0026  0.0009   224.68967   -972.41862  0:16
  20%  0.63 0.25 0.30 0.34 0.29 0.0037   10 111 P[5]=0.9885 0.0026  0.0009   210.93188   -972.34981  0:21
  25%  0.63 0.25 0.30 0.33 0.29 0.0038   10 111 P[5]=0.9892 0.0026  0.0009   224.49484   -972.24813  0:26
  30%  0.63 0.25 0.30 0.33 0.28 0.0035   10 111 P[5]=0.9901 0.0026  0.0009   216.40848   -972.27563  0:30
  35%  0.63 0.25 0.30 0.33 0.28 0.0033   10 111 P[5]=0.9907 0.0026  0.0009   224.00442   -972.50421  0:35
  40%  0.63 0.25 0.30 0.33 0.28 0.0031   10 111 P[5]=0.9912 0.0026  0.0009   213.49838   -972.43257  0:40
  45%  0.63 0.25 0.30 0.33 0.29 0.0034   10 111 P[5]=0.9905 0.0026  0.0009   209.30067   -972.37969  0:45
  50%  0.63 0.25 0.30 0.33 0.29 0.0033   10 111 P[5]=0.9906 0.0026  0.0009   232.27193   -972.51538  0:50
  55%  0.63 0.25 0.30 0.33 0.29 0.0035   10 111 P[5]=0.9900 0.0026  0.0009   203.35457   -972.50992  0:55
  60%  0.63 0.25 0.30 0.33 0.29 0.0037   10 111 P[5]=0.9895 0.0026  0.0009   206.99225   -972.50202  1:00
  65%  0.63 0.25 0.30 0.33 0.28 0.0038   10 111 P[5]=0.9898 0.0026  0.0009   217.40496   -972.48296  1:05
  70%  0.63 0.25 0.30 0.33 0.29 0.0038   10 111 P[5]=0.9893 0.0026  0.0009   216.82582   -972.50246  1:10
  75%  0.63 0.25 0.30 0.33 0.29 0.0037   10 111 P[5]=0.9896 0.0026  0.0009   204.99588   -972.56907  1:15
  80%  0.63 0.25 0.30 0.33 0.29 0.0036   10 111 P[5]=0.9897 0.0026  0.0009   217.51889   -972.54354  1:21
  85%  0.63 0.25 0.30 0.33 0.28 0.0037   10 111 P[5]=0.9893 0.0026  0.0009   212.80419   -972.51050  1:27
  90%  0.63 0.25 0.30 0.33 0.29 0.0037   10 111 P[5]=0.9896 0.0026  0.0009   223.34976   -972.47616  1:33
  95%  0.63 0.25 0.30 0.33 0.29 0.0036   10 111 P[5]=0.9893 0.0026  0.0009   215.63544   -972.48464  1:39
 100%  0.63 0.25 0.30 0.33 0.28 0.0037   10 111 P[5]=0.9895 0.0026  0.0009   201.02451   -972.47685  1:45

1:45 spent in MCMC

Summarizing the species-delimitation sample in file mcmc.txt

Number of species-delimitation models = 5

     model    prior    posterior
   1 000   0.200000   0.000000
   2 100   0.200000   0.000000
   3 101   0.200000   0.004630
   4 110   0.200000   0.005960
   5 111   0.200000   0.989410

Order of ancestral nodes:
  KCLH
  KC
  LH

Guide tree with posterior probability for presence of nodes:
((K, C)#0.995370, (L, H)#0.994040)#1.000000;;
```
The starting species delimitation model is generated by choosing at
random one of the models defined by the guide tree or starting tree.
After the chain has started, columns 2-6 give the
acceptance proportions for the conventional MCMC moves discussed above
while column 7 (rj) now shows the acceptance proportion for reversible jump proposals (collapsing 
or expanding species delimitations on the guide tree).
In general the larger this proportion, the more
efficient the rjMCMC algorithm is. However there is no optimal
acceptance proportion for the rjMCMC move, and a value close to 0 may
not necessarily mean a problem. If one model has posterior probability
close to 1, the acceptance proportion should be near 0 as well. Thus
both poor mixing of the rjMCMC algorithm and extreme posterior model
probabilities can cause the acceptance proportion for the rjMCMC to be
close to 0. It has been noted that if the rjMCMC algorithm is suffering
from poor mixing, different starting species trees often lead to
different results.
Next there are three numbers related to the rjMCMC move. 
The rest of the line shows the posterior means for $\theta$ and
$\tau_0$ of the root ancestral population (which exist in all
species-tree or species delimitation models), the current log MSC
density and average log sequence likelihood. The three numbers 
in columns 8-10 that are related to the rjMCMC move, for example,
```{ .yaml .no-copy }
10 111 P[5]=0.9895
```
represent the number of parameters in the current model (np), in this case 10, the species
delimitation, in this case 111, and the current posterior probability of the most probable
delimitation P[5]=0.9895. This means that delimitation model 5 (111) is the most favoured model,
with the posterior at 0.9895.

After the MCMC is finished, the program will summarize the sample. This
output is:
```{ .yaml .no-copy }
Summarizing the species-delimitation sample in file mcmc.txt

Number of species-delimitation models = 5

     model    prior    posterior
   1 000   0.200000   0.000000
   2 100   0.200000   0.000000
   3 101   0.200000   0.004630
   4 110   0.200000   0.005960
   5 111   0.200000   0.989410

Order of ancestral nodes:
  KCLH
  KC
  LH

Guide tree with posterior probability for presence of nodes:
((K, C)#0.995370, (L, H)#0.994040)#1.000000;;
```
The five delimitation models are listed, together with their posterior and prior probabilities.

#### MCMC output file
The MCMC output was saved in the file mcmc.txt. As the number of parameters changes when
the rjMCMC moves between models, the mcmc sample file may not be very
useful, so you can often ignore it. Right now the header line is generated
using the starting species tree and should be ignored. After the header
line, each line of output has the following numbers, separated by TABs:
iteration number, the number of parameters, the tree, the sampled
parameter values, and lnL. Here are the first 6 lines of mcmc.txt:
```{ .yaml .no-copy }
Gen	np	tree	theta_1K	theta_2C	theta_5KCLH	theta_6KC	theta_7LH	tau_5KCLH	tau_6KC	lnL
2	10	111	0.001207	0.003438	0.003433	0.000722	0.004404	0.003324	0.000834	0.000736	0.000516	0.000528	-972.035
4	10	111	0.001207	0.003438	0.001733	0.000874	0.004570	0.003914	0.001036	0.000736	0.000596	0.000528	-971.105
6	10	111	0.001207	0.003438	0.003411	0.001107	0.003353	0.001307	0.001944	0.000736	0.000708	0.000386	-973.526
8	10	111	0.001698	0.004529	0.001902	0.000913	0.003032	0.000349	0.001474	0.000607	0.000584	0.000346	-972.899
10	10	111	0.001979	0.005280	0.002217	0.001065	0.004159	0.000407	0.001990	0.000707	0.000680	0.000348	-970.674
```
If you know the unix command `grep`, you can retrieve the lines for the
same tree model to summarize the posterior for parameters in that model, for example,
```
    grep "3 111" mcmc.txt > result.Tree111.txt
```
In theory this should give you the same posterior as if you run analysis
A00 with the species tree fixed at tree 111. In practice it is
simpler to edit the Imap file and the control file to run analysis
A00.

**Running rjMCMC Algorithms.**

-   The rjMCMC algorithms for species delimitation allow the chain to
    move from one model to another but can have mixing problems. Check
    that similar results are produced from multiple runs using
    Algorithm0 and Algorithm1 regardless of the starting species tree.

    The starting tree is chosen by the program at random and will vary
    among runs. Make sure that some of your runs are started with the
    one-species model, some from the fully resolved tree, and some from
    other trees in between. If you get consistent results among runs
    with different starting trees and using the two algorithms, you are
    unlikely to have a convergence or mixing problem.

    If you have a computer with multicore, you can run those different
    combinations or replicates in different folders at the same time.

-   In medium-sized or large datasets with multiple loci, we have come
    across cases where the chain is stuck at the one-species model, or
    have difficulty moving into the one-species model. The problem is
    less common when there are only 1 or 2 loci. If you have a similar
    problem, you may start the analysis with 1 locus, then 2 loci, etc.
    to observe how the results change (you can do this by changing nloci
    in the control file as there is no need to change the sequence
    file).

-   The program used the burnin to adjust the step lengths for the
    proposals in the MCMC algorithm. If the rjMCMC stays in species
    model 0, which does not have any $\tau$ parameter, no information is
    collected during the burnin about the proposals to change $\tau$ and
    the automatically adjusted step length for $\tau$ can be very poor.

-   You probably need to evaluate the impact of the priors on $\theta$s
    and $\tau$. See the note about the gamma distribution later in this
    document.

### A11: Joint Species Delimitation and Species Tree Estimation
Analysis A11 `(speciedelimitation = 1, speciestree = 1)` assumes that the species delimitation is unknown
(e.g., population assignments as specified in the Imap file may not represent distinct species) and also that the species tree topology relating 
the species of any delimitation must be estimated (see [Yang and Rannala 2014](https://doi.org/10.1093/molbev/msu279)). The program calculates the joint posterior probabilities of different possible species delimitations and species tree topologies as well as of $\theta$s and $\tau$s. 

#### Input A11
Our example uses the dataset of frog sequence data generated by [Zhou et al 2012](https://doi.org/10.1111/j.1365-294X.2011.05411.x)
that we considered previously. The contents of the control file `A11.bpp.ctl` are shown below:
```
	seed =  -1

	seqfile = frogs.txt
	Imapfile = frogs.Imap.txt
	outfile = out.txt
	mcmcfile = mcmc.txt

	speciesdelimitation = 1 1 2 1  * species delimitation rjMCMC algorithm1 finetune (a m)
	speciestree = 1  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio
	speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 4  K  C  L  H
                    9  7 14  2
                 (((K, C), L), H);
         phase =   1  1  1  1
                  
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 5  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

   thetaprior = gamma 2 2000 # gamma(a, b) for theta (estimate theta)
     tauprior = gamma 2 1000 # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*     heredity = 1 4 4
*    locusrate = 1 5

      finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, 
mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 8000
      sampfreq = 2
       nsample = 100000
```
the lines:
```
	speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
	speciestree = 1  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio
```
specify that the rjMCMC algorithm 1, with $\alpha = 2$ and $m = 1$ in equations 6 and 7 of 
[Yang and Rannala 2010](https://doi.org/10.1073/pnas.0913022107) is used for species
delimitation and the species tree topology is being estimated
`speciestree = 1`. As noted above, the 3 parameters following the `speciestree` variable specify the probability that the
SNL move is used, rather than the subtree-pruning-regrafting (SPR) move to modify the species tree
and the expand/shrink ratios for the SNL move; these proportions can be adjusted on the interval $(0,1)$
and will affect mixing of the MCMC. The line:
```
	speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted
```
specifies the prior distribution on delimitations and tree topologies.
For analysis A11, `speciesmodelprior` can take the four values 0, 1, 2, 3, which
mean Priors 0, 1, 2, 3, respectively. As mentioned above, Prior 1 (which
is the default) assigns equal probabilities to the rooted species trees,
while Prior 0 means equal probabilities for the labeled histories
(rooted trees with the internal nodes ordered by age). Priors 2 and 3
assign equal probabilities for the numbers of species ($1/s$ each for 
$1,2,\ldots,s$ species given $s$ populations) and then divide up the
probability for any specific number of species among the compatible
models (of species delimitation and species phylogeny) either uniformly
(Prior 3) or in proportion to the labeled histories (Prior 2).
Priors 2 and 3 are discussed in [Yang and Rannala 2014](https://doi.org/10.1093/molbev/msu279).
Prior 3 may be suitable when there are many populations.

#### Output A11
When BPP is run using the above control file a summary of the progress of the MCMC is again printed to screen as follows:
```{ .yaml .no-copy }
     |            Acceptance proportions             |
Prgs | Gage Gspr thet  tau  mix   Sspr   Ssnl     rj | sp np        mldp mthet1   mtau1       log-PG        log-L
----------------------------------------------------------------------------------------------------------------
(some output omitted here)
   5%  0.65 0.26 0.30 0.28 0.30 0.0904 0.0000 0.0000    4 10 P[4]=1.0000 0.0041  0.0017   1239.61701  -4439.20352  0:59
  10%  0.65 0.26 0.30 0.27 0.30 0.0964 0.0000 0.0000    4 10 P[4]=1.0000 0.0042  0.0017   1221.22458  -4437.63784  1:35
  15%  0.65 0.26 0.30 0.26 0.30 0.0962 0.0000 0.0000    4 10 P[4]=1.0000 0.0042  0.0017   1238.23420  -4436.36392  2:13
  20%  0.65 0.26 0.30 0.26 0.30 0.0998 0.0000 0.0000    4 10 P[4]=1.0000 0.0042  0.0016   1207.06023  -4435.09525  2:54
  25%  0.65 0.26 0.30 0.26 0.30 0.1005 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1211.57636  -4434.31571  3:37
  30%  0.65 0.26 0.30 0.26 0.30 0.1018 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1199.97256  -4434.32883  4:19
  35%  0.65 0.26 0.30 0.26 0.30 0.1023 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1200.88102  -4433.94315  4:57
  40%  0.65 0.26 0.30 0.26 0.30 0.1002 0.0000 0.0000    4 10 P[4]=1.0000 0.0044  0.0016   1191.04407  -4434.06936  5:38
  45%  0.65 0.26 0.30 0.26 0.30 0.1045 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1193.71849  -4434.12002  6:21
  50%  0.65 0.26 0.30 0.26 0.30 0.1036 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1247.11644  -4434.33624  7:05
  55%  0.65 0.26 0.30 0.26 0.30 0.1017 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1218.84106  -4434.33556  7:43
  60%  0.65 0.26 0.30 0.26 0.30 0.0994 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1233.73293  -4434.00290  8:22
  65%  0.65 0.26 0.30 0.26 0.30 0.1005 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1235.85799  -4433.85179  9:06
  70%  0.65 0.26 0.30 0.26 0.30 0.0995 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1218.90538  -4434.20067  9:50
  75%  0.65 0.26 0.30 0.26 0.30 0.1002 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1214.40621  -4434.67870  10:32
  80%  0.65 0.26 0.30 0.26 0.30 0.0999 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1215.56021  -4435.19959  11:13
  85%  0.65 0.26 0.30 0.26 0.30 0.0980 0.0000 0.0000    4 10 P[4]=1.0000 0.0043  0.0016   1227.95558  -4435.47778  11:52
  90%  0.65 0.26 0.30 0.26 0.30 0.0963 0.0000 0.0000    4 10 P[4]=1.0000 0.0042  0.0016   1225.62882  -4435.73626  12:35
  95%  0.65 0.26 0.30 0.27 0.30 0.0958 0.0000 0.0000    4 10 P[4]=1.0000 0.0042  0.0017   1244.11759  -4435.87327  13:19
 100%  0.65 0.26 0.30 0.26 0.30 0.0949 0.0000 0.0000    4 10 P[4]=1.0000 0.0042  0.0017   1214.44552  -4436.19194  14:03

14:03 spent in MCMC


(A) List of best models (count postP #species SpeciesTree)
24746 0.247460 0.247460 4  (C H K L)  ((C, (H, L)), K);
20105 0.201050 0.448510 4  (C H K L)  (C, ((H, L), K));
10918 0.109180 0.557690 4  (C H K L)  (C, ((H, K), L));
 9358 0.093580 0.651270 4  (C H K L)  ((C, K), (H, L));
 6874 0.068740 0.720010 4  (C H K L)  ((C, (H, K)), L);
 5122 0.051220 0.771230 4  (C H K L)  (((C, H), L), K);
 4616 0.046160 0.817390 4  (C H K L)  (C, (H, (K, L)));
 3746 0.037460 0.854850 4  (C H K L)  (((C, K), H), L);
 3681 0.036810 0.891660 4  (C H K L)  (((C, H), K), L);
 3080 0.030800 0.922460 4  (C H K L)  (((C, L), H), K);
 1741 0.017410 0.939870 4  (C H K L)  ((C, L), (H, K));
 1693 0.016930 0.956800 4  (C H K L)  (((C, K), L), H);
 1691 0.016910 0.973710 4  (C H K L)  ((C, H), (K, L));
 1474 0.014740 0.988450 4  (C H K L)  ((C, (K, L)), H);
 1155 0.011550 1.000000 4  (C H K L)  (((C, L), K), H);

(B) 1 species delimitations & their posterior probabilities
100000 1.000000   4 (C H K L)

(C) 4 delimited species & their posterior probabilities
100000 1.000000 H
100000 1.000000 C
100000 1.000000 L
100000 1.000000 K

(D) Posterior probability for # of species
P[1] = 0.000000  prior[1] = 0.238095
P[2] = 0.000000  prior[2] = 0.238095
P[3] = 0.000000  prior[3] = 0.285714
P[4] = 1.000000  prior[4] = 0.238095
```
The arrangement of acceptance proportions, running means etc, are very similar to those of analysis
A10, except that there are two additional columns (with headers `Sspr` and `Ssnl`) that provide the
acceptance proportions for the SPR and SNL moves proposing changes to the species tree topology.
There are four sections of output that summarize the results.
Section (A) lists the best models in the decreasing order of the
posterior probabilities. Here a model is a full MSC model that species
both the species delimitation and species phylogeny. From this section,
one can easily construct the 95% or 99% credibility sets of models. This
section is further summarized to produce sections B, C, and D. Section
(B) gives the posterior probabilities for the top few delimitations.
Section (C) lists the delimited species and their posterior
probabilities, and section (D) lists the posterior probability for the
number of species, together with the prior probabilities calculated by
BPP.

#### MCMC output file
The sampled trees and delimitations from the MCMC are printed to the specified file `mcmc.txt`. An example of 5 lines from this file is: 
```{ .yaml .no-copy }
(C #0.010039: 0.001675, ((L #0.006075: 0.001466, H #0.003617: 0.001466) #0.002155: 0.000205, K #0.003031: 0.001671) #0.000955: 0.000004) #0.004713; 4
(C #0.009000: 0.001462, ((L #0.006482: 0.001228, H #0.001665: 0.001228) #0.001957: 0.000230, K #0.002632: 0.001458) #0.000511: 0.000004) #0.003204; 4
((C #0.009488: 0.001537, K #0.002525: 0.001537) #0.001082: 0.000004, (L #0.006926: 0.001384, H #0.002947: 0.001384) #0.001931: 0.000157) #0.003378; 4
((C #0.010464: 0.001806, K #0.002966: 0.001806) #0.001272: 0.000004, (L #0.008138: 0.001487, H #0.002791: 0.001487) #0.002480: 0.000323) #0.004458; 4
((C #0.009022: 0.001765, K #0.002899: 0.001765) #0.001243: 0.000004, (L #0.007955: 0.001449, H #0.003637: 0.001449) #0.002001: 0.000321) #0.004358; 4
```
Because the number of species delimited and the tree topology are both potentially changing during the MCMC this file is usually not useful for summarizing the results unless a custom script filtering the MCMC samples is used. This file is used by the program to produce the summaries printed to screen at the end of the run.

### Adjusting step lengths for MCMC moves (finetune)
The specified proposal step lengths are the initial values when automatic adjustment is
used and are the fixed step lengths when manual specification of step lengths is used.
The automatic adjustment of the finetune variable (MCMC step lengths) appears to be reliable
and is recommended. However, be sure to check that the resulting
acceptance proportions are neither too small nor too large. Here we provide
some guidance for manual adjustments of the step lengths. A good strategy is to use automatic
adjustments initially to generate good step lengths and then copy the step lengths into the
control file, while continuing to use automatic adjustment; that way the next time the control file is run the 
automatic adjustment procedure will start with near-optimal step lengths and less optimization
will therefore be required.

#### Adjusting the step lengths manually (finetune = 0)

First note the following line found in the control file `yu2001.bpp.ctl`.
Here finetune = 1 means that MCMC step lengths will be adjusted
automatically, and the specified values are used as initial values.
```
    finetune = 1: 2 0.00002 0.0001 0.0005 0.5 0.2 1.0
```
There are seven finetune steplengths here. They are in a fixed order and
always read by the program even if the concerned proposal is not used.
The first five of them are $\epsilon_1$, $\epsilon_2$, $\epsilon_3$,
$\epsilon_4$, and $\epsilon_5$, described in 
([Rannala and Yang 2003](https://doi.org/10.1093/genetics/164.4.1645)). These are
the step lengths used in the MCMC proposals that (1) change internal
node ages in the gene tree, (2) prune and re-graft nodes in the gene
tree (SPR), (3) update $\theta$, (4) update $\tau$s using the rubber-band
algorithm, and (5) implements the mixing step. The 6th and 7th are for
the proposals that change the locus rates or heredity multipliers and
that change the sequencing errors, respectively. If the model assumes
the same rate for all loci and does not use heredity multipliers, the
6th proposal step is not used. If the model assumes no sequencing
errors, the 7th step is not used. The acceptance proportions for the
first five proposals are always printed out on the screen, but those for
the 6th and 7th are printed out only if the concerned proposal is used
in the model. In the example above, only the first five proposals are
used in the algorithm and we will manually change the step lengths so that the
acceptance proportions become close to 30%. 
First, edit the line, changing it to:
```
    finetune = 0: 2 0.00002 0.0001 0.0005 0.5 0.2 1.0
```
The basic strategy for adjusting fine-tune parameters is if the acceptance proportion
is too small (e.g., $<0.10$), decrease the corresponding finetune
parameter. If the acceptance proportion is too large (e.g., $>0.80$),
increase the finetune parameter.
Run the program for a small number of iterations and look at the screen
output to examine the acceptance proportions:
```
   0%  0.71 0.19 0.34 0.00 0.41   0.0004    463.81031  -12720.98567  0:04
   5%  0.72 0.19 0.32 0.00 0.41   0.0003    458.65988  -12720.53426  0:04
  10%  0.72 0.19 0.34 0.00 0.42   0.0004    454.99289  -12720.66361  0:05
  15%  0.72 0.19 0.34 0.00 0.41   0.0004    481.20520  -12720.78839  0:06
```
Here the first acceptance proportion, at 0.72, appears too large, which
means that the corresponding finetune parameter (2 in the control file) is
too small, similarly the second acceptance proportion, at 0.19, appears 
too small, which means that the corresponding finetune parameter is too large. 
Terminate the run (Ctrl-C) and edit the finetune values in the
control file and then run the program again, checking the acceptance rates. 
Repeat this process a few times until every
acceptance proportion is neither too small nor too large. In this
example, changing the second finetune parameter (SPR) from $0.00001$ to 
$0.000003$ brings the acceptance proportion to about $0.30$, namely:
```
  10%  0.71 0.31 0.34 0.00 0.30   0.0004    447.59035  -12721.04474  0:06
  15%  0.71 0.30 0.33 0.00 0.30   0.0004    487.91342  -12720.96697  0:07
  20%  0.71 0.30 0.32 0.00 0.30   0.0003    430.67545  -12720.83600  0:08
```
However, increasing the first finetune parameter (even by an order of magnitude to 20) has little effect. Sometimes it will not
be possible to reduce the acceptance proportion for a parameter, as in
this case (even using automatic adjustment produces an acceptance proportion of about $0.71$ in this case).

The MCMC proposals discussed above are used in all four analyses (A00, A01, A10, A11),
so that the description here applies to all of them. Note that the
finetune parameters affect the efficiency of the MCMC or how fast one
can obtain reliable results. In theory they do not change the results if
runs using different finetune parameters (with potentially different acceptance proportions) are all run for long enough to generate
reliable results.


## Simulating Gene Trees and Alignments Using BPP

The `--simulate` option of BPP4 replaces the earlier MCCOAL program distributed with
BPP versions 3.4 and earlier. This option can be used to simulate gene trees and
sequence alignments at multiple loci under the MSC ([Rannala and Yang 2003](https://doi.org/10.1093/genetics/164.4.1645)) and
MSC-I ([Flouri et al 2020](https://doi.org/10.1093/molbev/msz296)) models.
The simulation program allows a GTR+G DNA substitution model (and its special cases) 
and allows among-loci heterogeneity in the process of sequence evolution. The different loci can have
different exchangeability parameters ($a, b, c, d, e, f$) and base
compositions in the GTR model, different overall evolutionary rates,
differing extents of among-site rate variation (reflected in the alpha
parameter for the gamma model), and different among-species rate drift processes
(violations of the clock). See @Shi2018 for examples of such simulations.
The simulation option also includes a MSC-M continuous migration model, with a
user-specified migration rate matrix between species/populations.

### Running the simulation program
To run BPP using the simulation option, create a suitable control file (for example, named MCcoal.ctl) and type the following:
```
bpp --simulate MCcoal.ctl
```
Always examine the screen output carefully to confirm that the program has read and executed the control file correctly.

### Simulation under the MSC model
#### The control file
Here we consider a simple example of a simulation under the MSC model with no migration or introgression.
The example control file for simulating sequence alignments and gene trees has the following content:
```
    seed = 12345
    seqfile = MySeq.txt  * comment out this line if you don't want seqs
    treefile = MyTree.tre  * comment out this line if you don't want trees
    Imapfile = MyImap.txt
    *   concatfile = concatenatedfile.txt  * concatenated alignment
    modelparafile = modelparas.txt   * comment out this line if you don't want seqs

    species&tree = 4  A  B  C  D
    3  2  1  1
    ((A #0.01, B #0.01) :0.01 #.01, (C, D) :0.011 #0.01) : 0.012 #0.01;
    phase =    1  1  0  0   * 0: do not phase (fully phased seqs), 1: diploid unphased seqs

    loci&length = 100 1000 * number of loci & number of sites at each locus
    *    locusrate = 0 (default: same rate for all loci)
    *    locusrate = mu_bar a_mui prior
    *    clock = 1 (default: strict clock)
    *    clock = 2 v_bar a_vi prior dist  (independent-rates)
    *    clock = 3 v_bar a_vi prior dist  (correlated-rates)

    model = 7       * model: 0:JC69, 7:REV (GTR)
    Qrates = 0  10 5 5 5 5 10     * 1: fixed; 0: dirichlet, for TC TA TG CA CG AG
    basefreqs = 0  10  10  10  10   * 0: random, Dirichlet(aT,aC,aA,aG), for base frequencies
    alpha_siterate = 0  100  20  5  * G(a, b) for alpha for sites & K for discrete gamma
```
We now go through each control file variable and explain its purpose.

*seed*  
Use `seed = -1` to simulate different datasets each time the program is
run (obtaining a random number seed from the computer clock). 
If you use a positive integer as a seed the data will be exactly the same
each time the program is run with that control file; this is usually not
what you want to do.

*seqfile, treefile, Imapfile, concatfile, modelparafile*  
These variables specify the names of files to be generated. If you want to simulate gene trees and
not sequence files, you can comment out the line for the `seqfile`.
The file specified by `concatfile` will contain the sequence alignment concatenated across loci.
The sequence files can become very large. The file specified by `modelparafile` contains the
parameter values for the GTR model, as well as the locus rate and the
rates for nodes on the species tree. There is one line of output per
locus.

*phase*  
This variable is used in the same way as in the BPP inference (data
analysis) mode. Each species is assigned a switch/flag, with 0 meaning
fully phased (haplotype) sequences and 1 for unphased diploid sequences.
With the specifications in the example (sampling configuration 3 2 1 1
and phase = 1 1 0 0), the program will simulate 6, 4, 1, 1 sequences for
A, B, C, D, respectively, and then combine every two sequences from A,
B, C into one single diploid sequence with heterozygote sites
represented using ambiguity characters (YRMKSW). If you use the same
random number seed, the diploid data generated using the sampling
configuration (3 2 1 1) with phase = 1 1 0 0 should match the fully
resolved sequences generated using the sampling configuration (6 4 1 1)
with phase = 0 0 0 0. The program also prints out the fully resolved
sequences in a file with the suffix `_full` (or `MySeq_full.txt` in the
example).

*loci&length*  
This variable specifies the number of loci to be simulated and the length 
(number of sites) for each loci. In this example, 100 loci,
each of 1000 sites, will be simulated.

*species&tree*  
This block specifies the number and names of species,
the species tree, and the parameters under the multispecies coalescent
model, including the species divergence times ($\tau$s) and population
size parameters ($\theta = 4N\mu$). For the example above, 100 loci,
each of 1000 sites, will be simulated, on the species tree ((A, B), (C, D)), 
with 3, 2, 1, 1 sequences for A, B, C, D, so that there are 7
sequences at each locus. The divergence time parameters ($\tau$s) are
after `:` in the tree, while the population size parameters ($\theta$s)
are after `#`. Thus we have 0.01 for all $\theta$s, and
$\tau_{ABCD} = 0.012$, $\tau_{AB} = 0.01$, and $\tau_{CD} = 0.011$. We
need $\theta_A$ and $\theta_B$ because 2 or more sequences are sampled
from species A and B. Parameters $\theta_C$ and $\theta_D$ are
unnecessary in this case as only one sequence is sampled from C and D:
if you specify them, they will be ignored by the program. Note that both
$\theta$s and $\tau$s are measured by the expected number of mutations
per site, so that $\theta = 0.01$ means that two sequences sampled from
the population are 1% different. For reference, the estimate for the
human species is $\theta_H = 0.0006$.

*model*  
This variable can take two possible values: 0 for JC69 and 7 for GTR.
The next few lines are relevant for the GTR model.

*Qrates*  
This variable specifies the exchangeability parameters in the GTR model
($a, b, c, d, e, f$ for TC, TA, TG, CA, CG, and AG, in Yang 1994).
There are two options, either to have those rates fixed for all loci or
to have them sampled at random. The first has the following format (note
that the first 0-1 number is a switch):

    Qrates = 1  2 1 1 1 1 2   *  1: fixed; for TC TA TG CA CG AG

The input is one integer value (0 or 1) which acts as a switch followed
by six real numbers. The line above fixes the rates at
$a = 2, b = c = d = e = 1$ and $f = 2$, so that the transition rate is
twice as high as the transversion rate, and the model corresponds to K80
or HKY. Note that only the relative rates matter, because the rate
matrix ($Q$) is scaled so that the average rate (over the base
frequencies) is 1 and branch lengths are measured in the expected number
of mutations/substitutions per site. Nevertheless the program expects
six rate parameters in the input here:

    Qrates = 0  10 5 5 5 5 10   *  0: dirichlet, for TC TA TG CA CG AG

The second option, above, is to sample $a, b, c, d, e$ from the
Dirichlet distribution for every locus. The above specifies Dir(10, 5,
5, 5, 5, 10), with parameters
$\alpha_a = 10, \alpha_b = \alpha_c = \alpha_d = \alpha_e = 10$, and
$\alpha_f = 10$. The rates generated from the Dirichlet sum to 1, and
they are rescaled. Note that larger values for those parameters mean
less variance, while the mean for $a$, say, is given by
$\alpha_a/\alpha$ with
$\alpha = \alpha_a + \alpha_b + \alpha_c + \alpha_d + \alpha_e +
\alpha_f$. The specification here gives on average a
transition/transversion rate ratio of 2 (if the base frequencies are
fixed at $\frac{1}{4}$ each), but it varies among loci.

*basefreqs*  
The base frequency parameters
($\pi_T, \pi_C, \pi_A, \pi_G$) can also be fixed or sampled from the
Dirichlet distribution, like the `Qrates`. The two options are as
follows:
```
    basefreqs = 1 0.15 0.35 0.15 0.35  * 1: fixed; Base frequencies are for TCAG.
    basefreqs = 0  10  10 10  10       * 0: random, Dirichlet(aT, aC, aA, aG), for base frequencies
```

*alpha_siterate*  
The gamma shape parameter ($\alpha$) for rate
variation among sites of the same locus can be fixed or sampled from a
gamma distribution. The possible options are as follows. When alpha is
fixed, the same value is used for all loci. Otherwise different
$\alpha$s are used for different loci. Note that given alpha, the rates
for sites at the same locus have a G(($\alpha, \alpha$) distribution
with mean 1 [@Yang1993; @Yang1994]. The fourth option below will allow
the program to sample $\alpha$ from $G(100, 20)$, with mean 5, for each
locus and then use $G(\alpha, \alpha)$ to describe the among-site rate
variation for the locus:
```
    alpha_siterate = 1 0            * 1: alpha fixed at 0(inf), one rate for all sites at the locus
    alpha_siterate = 1  5.6 0       * 1: alpha fixed at 5.6, K = 0(inf) for continuous gamma
    alpha_siterate = 1  5.6 5       * 1: alpha fixed at 5.6, K = 5 categories for discrete gamma
    alpha_siterate = 0  100  20  5  * 0: alpha sampled from G(100, 20), with mean 5, K = 5 for discrete gamma
```

The variables `locusrate` and `clock` are used to specify
relaxed-clock models, with rate variation among branches and among loci.
These are the same models as implemented in the inference program.
Please see [Models of Substitution Rate Variation](#models-of-substitution-rate-variation) for details.

*locusrate*  
This specifies variable rates for different loci. The default is
0, meaning that all loci have the same rate. To specify variable rates
for loci, use the following format:
```
    locusrate = 0 (default)
    locusrate = mu_bar a_mui prior
    locusrate = mu_bar a_mui
    locusrate = 1.0 5.0 dir
```
Here `mu_bar` is the mean rate across loci.
If `prior = dir` (gamma-Dirichlet), the total rate ($L\bar\mu$) is
partitioned to generate the locus-rates ($\mu_i$) according to the
Dirichlet distribution with concentration parameter $\alpha$ (`a_mui`).
Large $\alpha$ (10 or 100, say) means that the rates are similar among
loci, while small values (1 or 0.5, say) mean that rates are highly
variable among loci. If `prior = iid` (conditional i.i.d.), the locus
rates ($\mu_i$) are i.i.d. given the mean rate ($\bar\mu$):
```
    mu_i ~ Gamma(a_mui, a_mui/mubar)
```

*clock*  
This variable is used to specify strict-clock or relaxed-clock models.
The default is `clock = 1` for strict clock, while `clock = 2` is the
independents-rates model, and `clock = 3` is the correlated-rates model
(which is implemented for the MSC model with no introgression and is
unavailable for the MSC-I model).
```
    clock = 1 (default: strict clock)
    clock = 2 v_bar a_vi prior dist  * independent-rates model
    clock = 2 0.1 10 dir G    * gamma-dirichlet for loci & gamma for branches
    clock = 2 0.1 10 iid LN   * iid for loci & log-normal for branches
    clock = 3 v_bar a_vi prior dist  * correlated-rates model
    clock = 3 0.1 10 iid G    * iid for loci & gamma for branches
    clock = 3 0.1 10 iid LN   * iid for loci & log-normal for branches
```
The specification (`clock = 2 0.1 10 dir G`) means the following. First
the average variance parameter is = 0.1. Then the total $L\bar\mu$ is
partitioned into $\mu_i$ (for $i = 1, 2, \cdots, L$) according to the
Dirichlet distribution with concentration parameter $\alpha = 10$
(`a_vi`). In contrast, the specification (`clock = 2 0.1 10 iid LN`) means that
the variance parameter for locus $\nu_i$ is generated as
$\nu_i|\bar\nu \sim G(\alpha, \alpha/\bar\nu)$ with shape parameter
$\alpha = 10$ (`a_vi`) and mean $\bar\nu = 0.1$.
In both prior models (`prior=dir` and `iid`), $\alpha$ (`a_vi`) is
inversely related to the variance in $\mu_i$ among loci: use small
values of $\alpha$ (2, 1, or 0.5) if the clock nearly holds at some loci
but is seriously violated at others, and large values (such as 10 or
100) if the clock is violated in the same extent at different loci.

`clock = 2`  
Given the locus rate $\mu_i$ (specified using the
`locusrate` variable) and the variance parameter $\nu_i$ (specified
using the `clock` variable) for locus $i$, the rates for species-tree
branches are specified as follows. If `dist = G` (for gamma), the rate
for branch $j$ at locus $i$ has the following gamma distribution with
mean $\mu_i$ and variance $\nu_i$

$r_{ij} | \mu_i, \nu_i \sim G(\mu_i^2/\nu_i, \mu_i/\nu_i).$

If `dist = LN` (for log-normal), the rate has the log-normal
distribution 

$r_{ij} | \mu_i, \nu_i \sim LN(\mu_i, \nu_i).$

`clock = 3`  
The correlated-rates model is specified similarly, and for
the MSC model only. Given the locus rate $\mu_i$ (specified using the
`locusrate` variable) and the variance parameter $\nu_i$ (specified
using the `clock` variable) for locus $i$, the rates for species-tree
branches are specified recursively, starting from the root towards the
tips.
If `dist = G` (for gamma), the rate for branch $j$ at locus $i$ has the
following gamma distribution with shape parameter
$\alpha = \mu_i^2/\nu_i$ and with the mean to be the rate for the
parental branch of $j$ or $r_{i,p(j)}$:

$r_{ij} | \mu_i, \nu_i, r_{i,p(j)} \sim G(\frac{\mu_i^2}{\nu_i}, \frac{\mu_i^2}{\nu_i r_{i,p(j)}}).$

If `dist = LN` (for log-normal), the rate for branch $j$ at locus $i$
has the log-normal distribution with the mean to be the parental rate
and variance parameter $\nu_i$:

$r_{ij} | \nu_i, r_{i,p(j)} \sim LN(r_{i,p(j)}, \nu_i).$

After the rates for branches (populations) on the species tree are
generated for the locus, they are used to calculate the gene-tree branch
lengths. Gene-tree branches residing in the same species/populations
have the same rate (the rate of that species at the locus), while those
residing in different species may have different rates. A branch on a
gene tree may have segments residing in different species/populations,
and the length of the branch (measured by the expected number of
mutations/substitutions per site) is calculated by adding up those
segments.

#### Program output
Running the program on the above control file produces the following screen output, indicating
that the simulation has been successful:
```{ .yaml .no-copy }
bpp --simulate MCcoal.ctl 
Site rates: alpha sampled from Gamma(100.000000,20.000000), K = 5
bpp v4.6.2_linux_x86_64, 15GB RAM, 8 cores
https://github.com/bpp/bpp

Detected CPU features: mmx sse sse2 sse3 ssse3 sse4.1 sse4.2 popcnt avx avx2
Auto-selected SIMD ISA: AVX2

Substitution model: GTR with generated rates from Dirichlet
                    10.000000 5.000000 5.000000 5.000000 5.000000 10.000000
Base frequencies: generated from Dirichlet
                  a_T=10.000000 a_C=10.000000 a_A=10.000000 a_G=10.000000
4 species: A (3) B (2) C (1) D (1)
  ((A #0.010000: 0.010000, B #0.010000: 0.010000) #0.010000: 0.002000, (C: 0.011000, D: 0.011000) #0.010000: 0.001000) #0.010000;

Map of populations and ancestors (1 in map indicates ancestor):
  Species   1  2  3  4  5  6  7
1 A         1  0  0  0  1  1  0  tau = 0.000000  theta = 0.010000
2 B         0  1  0  0  1  1  0  tau = 0.000000  theta = 0.010000
3 C         0  0  1  0  1  0  1  tau = 0.000000  theta = 0.000000
4 D         0  0  0  1  1  0  1  tau = 0.000000  theta = 0.000000
5 ABCD      0  0  0  0  1  0  0  tau = 0.012000  theta = 0.010000
6 AB        0  0  0  0  1  1  0  tau = 0.010000  theta = 0.010000
7 CD        0  0  0  0  1  0  1  tau = 0.011000  theta = 0.010000

Sequence data file -> MySeq.txt
Trees -> MyTree.tre
Model parameters for loci -> modelparas.txt
Tags to species mapping (Imap) -> MyImap.txt
```
Several output files are produced. The file `MySeq.txt` contains the simulated sequence alignments in Phylip/BPP format:
```{ .yaml .no-copy }
7 1000 

a1^A        ACATGCCACT ACACCGGACA TTAGAGCTAC GCGAGCTCTT TGWACTCGTT ACAACCTCAC AATGTAATAG TGTAATTACT TAGCTTAAGT AMT
CTTGATC TTTTCTATAR TTACATCGAA CGAATATCTG TACCTGAGCT TTTCACGATC CCTYCCTTCA TCTAGCAATC GTAAACTWTG TAAATTTCCT ACAGTCA
TAY TRATAACCCA TATCATGTAA CCAGACCCCC CTWAATTTMT CGAGAATGCC CMTCTTCATT CAATCCCAAA CAGTCAAAGC GTCCACAAAT TAATTGCACA 
TAACCATCGT AATTGGAGAA ATTCCACCAT TCAATAGTTT ACTTGATTMA AACTAGTGGA TCCCATGTAA GCATAATCCC CCTACYCCAT ATGTATCTTA GAGT
MTCGTT TGRGACTTTM CATGCCCCAG CTGTTTTTAA ATGCCAGTAA ACARCGGCGA CTACTATGTC TTACGACAGA ACGAGGTTTC CACTGAGAGT RCAAAGGA
TC ACTCGTCTAT AGTTCAAAAW TCATTGCACG ATGTTCTCGA TCTACCTAAA TCTTGCCTGC TTTTTCMCCA TGCAATACAA CACCGTTAAA CAATATCTTA A
TATCAATCG TACAGTCGAT GCAGTATAGC CGATACTCCT TAACCTAACT YTTCAATATG ACMTAATACT CAACTGTTGC GRTACCAATT ATCCMTGCTG TACTG
ATAAC TTTATTGAAT GCCAATMCTA TWATGATCTA TTTCAGAYTT AATTCACAAT GAATTCAATA CGAATTTTCC TATCGTCAGS ACATCAATTA GATTCACAT
G TKGRAAGATC TATTTTGCAT CGACAAATGC ACCACCGCTC CCTCGCCTCT AAATCGCTCA CAATCGAATT AGAAAATATA ACACWTCCAT TGATCAASCT AA
AAGCGCAC AGGGAGCGAT TGGTTCGAAA ATGAGTTACA ACGTAAAGTT CTCAGCATAA TGGRGGCTAT ACACCCAAAA
a2^A        ACATGCCACT ACACCGGACA TTAGAGCTAC GCGAGCTCTT TGTACTCGTT ACAACCTCAC AATGTAATAG TGTAATTACT TAGCTTAAGT AAT
CTTGATC TTTTCTATAA TTACATCGAA CGAATATCTG TACCTGAGCT TTTCACGATC CCTCCCTTCA TCTAGCAATC GTAAACTTTG TAAATTTCCT ACAGTCA
TAT TAATAACCCA TATCATGTAA CCAGACCCCC CTTAATTTCT CGAGAATGCC CATCTTCATT CAATCCCAAA CAGTCAAAGC GTCCACAAAT TAATTGCACA 
TAACCATCGT AATTGGAGAA ATTCCACCAT TCAATAGTTT ACTTGATTCA AACTAGTGGA TCCCATGTAA GCATAATCCC CCTACCCCAT ATGTATCTTA GAGT
CTCGTT TGGGACTTTC CATGCCCCAG CTGTTTTTAA ATGCCAGTAA ACAGCGGCGA CTACTATGTC TTACGACAGA ACGAGGTTTC CACTGAGAGT GCAAAGGA
TC ACTCGTCTAT AGTTCAAAAT TCATTGCACG ATGTTCTCGA TCTACCTAAA TCTTGCCTGC TTTTTCCCCA TGCAATACAA CACCGTTAAA CAATATCTTA A
TATCAATCG TACAGTCGAT GCAGTATAGC CGATACTCCT TAACCTAACT CTTCAATATG ACATAATACT CAACTGTTGC GATACCAATT ATCCATGCTG TACTG
ATAAC TTTATTGAAT GCCAATCCTA TTATGATCTA TTTCAGATYT AATTCACAAT GAATTCAATA CGAATTTTCC TATCGTCAGC ACATCAATTA GATTCACAT
G TGGGAAGATC TATTTTGCAT CGACAAATGC ACCACCGCTC CCTCGCCTCT AAATCGCTCA CAATCGAATT AGAAAATATA ACACTTCCAY TGATCAAGCT AA
AAGCGCAC AGGGAGCGAT TGGTTCGAAA ATGAGTTACA ACGTAAAGTT CTCAGCATAA TGGGGGCTAT ACACCCAAAA
...
```
The file MyTree.tre contains the simulated gene tree for each locus in Newick format:
```{ .yaml .no-copy }
((((((a1a^A:0.000217,(a2b^A:0.000182,a3a^A:0.000182):0.000035):0.000362,a2a^A:0.000579):0.000065,a3b^A:0.000644):0
.011379,c1^C:0.012023):0.001408,(a1b^A:0.010443,(((b1a^B:0.001294,b2b^B:0.001294):0.000858,b1b^B:0.002152):0.00144
4,b2a^B:0.003596):0.006847):0.002989):0.010039,d1^D:0.023471):0.000000; [TH=0.023471, TL=0.091502]
(((((a1a^A:0.000584,((a2b^A:0.000357,a3b^A:0.000357):0.000057,a3a^A:0.000414):0.000170):0.001213,a1b^A:0.001797):0
.000482,a2a^A:0.002279):0.008361,(b1a^B:0.001291,(b1b^B:0.000558,(b2a^B:0.000345,b2b^B:0.000345):0.000213):0.00073
2):0.009350):0.010901,(d1^D:0.014926,c1^C:0.014926):0.006615):0.000000; [TH=0.021541, TL=0.076273]
((((a1a^A:0.001025,((a1b^A:0.000262,a2a^A:0.000262):0.000003,a2b^A:0.000265):0.000760):0.003850,(a3a^A:0.003915,a3
b^A:0.003915):0.000961):0.005929,((b1a^B:0.000958,(b2a^B:0.000139,b2b^B:0.000139):0.000819):0.000315,b1b^B:0.00127
3):0.009531):0.005512,(d1^D:0.015947,c1^C:0.015947):0.000369):0.000000; [TH=0.016316, TL=0.072093]
...
```
The file `MyImap.txt` contains the Imap data. These files can be used as input to inference programs such as
BPP along with a suitable control file that matches the simulated data
.
### Simulation with migration under the MSC-M model
Here we consider a simple example of a simulation under the MSC-M model with continuous migration.
#### The control file
The example control file for simulating sequence alignments and gene trees has the following content:
```
seed = -1
treefile = mytree.tre
Imapfile = myimap.txt
seqfile = mydata.txt
*concatfile = concat.txt
modelparafile = modelparas.txt
#phase = 0 0 0
species&tree = 3
A B C
4 4 4
((A #0.015, B #0.025) S #0.015 :0.01, C #0.025) :0.02 #0.025;
loci&length = 2000 500
model = 0
* JC
migration = 2
S C 0.1
A C 0.2 5
```
There are two major differences from the previous (MSC) example. First, in the species tree
the ancestral population to clade (A, B) is given a label S:
```
species&tree = 3
A B C
4 4 4
((A #0.015, B #0.025) S #0.015 :0.01, C #0.025) :0.02 #0.025;
```
A label is needed if migration is going to be specified for the ancestral population.
Second, there is a new `migration` block:
```
migration = 2
S C 0.1
A C 0.2 5
```
Here, two migration connections are specified with `migration = 2`. The first migration connection is from source population 
S (ancestor of AB) to target population C and migration occurs at the rate $M_{SC} = 0.1$ migrants per generation. The second 
migration connection is from source population A to target population C. In this case,
the migration rate $M_{AC}$ varies among loci with the mean to be 0.2 and the shape parameter to be 5.
The migration rate from populations $i$ to $j$ is defined as $M_{ij} = N_j m_{ij}$ , where $m_{ij}$ is the
proportion of individuals in population $j$ that are migrants from population $i$ in each generation. In the example above,
$M_{SC} = 0.1$, which means that on average 0.1 individuals are immigrants from population S
to population C. We use forward-time for our definition of rates as the rates are then more biologically meaningful.
This control file format is used in BPP version 4.4 onward, and is different from the old MCcoal format, in which the labels and 
order for populations are fixed by the program. The control file included in the release ( `MCcoal.im-3s-saturated.ctl` ) simulates data on a
species tree for three species, with eight migration rates.


## Advanced Features of BPP

### Threading and Checkpointing


#### Threads

BPP4 can use Single-Instruction-Multiple-Data (SIMD)
instruction sets (also called vector instructions) available on modern
processors. It automatically detects and uses the best instruction set
available so you don't have to do anything about this. However, you can
manually specify and force an instruction set in the control file using
the arch tag. For instance, arch=SSE forces BPP to use the
SSE instruction set even if AVX (which has twice the register size of
SSE and thus in theory yields twice the speedup) is present on the
system. To completely disable vector instruction, one can use the
arch=CPU option. The available values for the arch tag are 'CPU', 'SSE',
'AVX' and 'AVX2'.

![Fig-Numa](https://github.com/bpp/bpp-manual/blob/191d107413566b145360597fdfdfef3e68ce9173/figures/fig-numa.png?raw=true)

Modern computer architecture using Non-Uniform momery access (NUMA).

BPP4 can use multiple cores/threads. The control variable threads has the following syntax.
```
    Threads = 8 19 1
```
This means BPP will use 8 threads, starting from
core/thread 19, with increment 1; in other words hardware threads 19-36
will be used. The software threads are pinned onto the hardware threads
and they not allowed to migrate. We found that this helps with
performance.

Modern clusters based on microprocessors use a common shared-memory
configuration called NUMA (for non-uniform memory access). The cluster
typically consists of two or four microprocessors interconnected on a
local bus to a shared memory on a single motherboard.

You can use commands like `lscpu` to see the machine architecture. For
example, on our server lscpu gives the following output:
```
    Architecture:          x86_64
    CPU op-mode(s):        32-bit, 64-bit
    Byte Order:            Little Endian
    CPU(s):                144
    On-line CPU(s) list:   0-143
    Thread(s) per core:    2
    Core(s) per socket:    18
    Socket(s):             4
    NUMA node(s):          4
    ...
    NUMA node0 CPU(s):     0-17,72-89
    NUMA node1 CPU(s):     18-35,90-107
    NUMA node2 CPU(s):     36-53,108-125
    NUMA node3 CPU(s):     54-71,126-143
```
This means that hardware threads 1-18 are on CPU0, while threads 73-90
are the hyperthreads on CPU0. The specification Threads = 8 19 1 thus
uses 8 threads on CPU1. Our recommendation is that you do some tests
yourself using threads = 2 or 4 or 8, and use options that work well for
your servers. Adjust the step lengths first. Then use burnin=0 and a
small nsample so that the running time is about 1-2 minutes. Change
threads and record running time.

Make sure all threads for the same BPP job are on the same
CPU. Use lscpu to see the machine configuration and htop to examine the
load on the hardware threads.

#### Checkpointing

This option instructs BPP to create a checkpoint file
(save the progress) after a specified number of MCMC iterations. The
MCMC loop in BPP consists of
`N = burnin + sampfreq*nsample` iterations, numbered as $1,2,\cdots,N$.
There are two formats, where X, Y are whole numbers:
```
    checkpoint = X
    checkpoint = X Y
```
In the first format, a single checkpoint is created after X MCMC
iterations (including burnin).

In the second format, a checkpoint is created after X steps, and then
additional checkpoints are created every Y steps. For example, if X =
10000 and Y = 10000, the first checkpoint is created at MCMC iteration
10000, and further checkpoints at iterations 20000, 30000, etc.

The checkpoint files are named OUTFILE.Z.chk, where OUTFILE is the label
specified for the outfile option, and Z is the number of the checkpoint
file, starting from 1 and incrementing with every new checkpoint file.

To resume from a checkpoint file, use the resume switch, i.e.
```
    bpp --resume checkpoint-file.chk
```

### Marginal Likelihood Calculations

The C program BFdriver generates control files and job subscription
scripts for running MCMC (using BPP or MCMCtree) to calculate the
marginal likelihood (or the Bayes factor), as described in [Rannala and Yang 2017](https://doi.org/10.1093/sysbio/syw119).

The program takes a control file you provide (such as bpp.ctl) and
generates $K = 16$ control files with different beta values, which are
used to run bpp to sample from the different power posterior
distributions. The program also generates job submission scripts and
submit the jobs using qsub. All generated control files and output files
are in the same current directory. The frogs dataset in the bpp release
(([Yang 2015](https://doi.org/10.1093/czoolo/61.5.854)) is used as the example.

You need the following: a linux system with SUN grid engine managing job
submission (including commonds such as `qsub`, `qstat`, `qdel`, etc.),
and a C compiler. If you don't have this job submission system, you can
use BFdriver to generate the control files and run the MCMC jobs from
the command line.

#### Compiling and running BFdriver
```
    cc -o BFdriver -O3 BFdriver.c tools.c -lm
    BFdriver <controlfilename> <npoints> <scriptname.sh>
    BFdriver A00.ctl 16 tmp.sh
```
You may need to edit the following two lines inside BFdriver.c, and if
you do, remember to compile the program after editing. Here bpp is
assumed to be on your search path. You can use a full path for the
executable prorgam, such as `/home/gooduser/bin/bpp3.3`. Also the second
line is for submitting the jobs using `qsub`. Here the limits are set to
4G of RAM and 360 hours of running time. Check those values if necessary
(and recompile).
```
    fprintf(fcommand, "     echo \"bpp %s.b$I.ctl > log.b$I.txt\" > %s\n", ctlf, scriptf);
    fprintf(fcommand, "     qsub -S /bin/bash -l h_vmem=4G -l tmem=4G -l h_rt=360:0:0 -cwd %s\n", scriptf);
```

#### Running the program 

Create a folder inside `frogs/`, say `bf1`:
```
    mkdir frogs/bf1
    cd frogs/bf1
```
Prepare a control file (A00.ctl, say) for the A00 analysis in the
folder. Check that it works. This should have a fixed species tree,
which is ((K, C), (L, H)). Species delimitation and species tree
estimation should be turned off. The control file specifies the priors
for theta, tau, and also specifies burnin, nsample, sampfreq etc. Run
bpp to confirm that the control file works. Then run BDdriver as
follows:
```
    BFdriver A00.ctl 16 tree1.sh
```
Here A00.ctl is the control file we have prepared. $K = 16$ is the
number of points in the Gauss-Legendre qradrature algorithm for
numerical integration. You can use 8 for testing, and 16 or 32 for real
calculation. `tree1.sh` is the temporary script file for job submission
using `qsub`. The BFdriver command does a few things. First it reads the
control file specified (`A00.ctl`) and creates 16 control files with
names like `A00.b01.ctl`, \..., `A00.b16.ctl`. Each of those control
files has the same content as `A00.ctl` except that one extra line is
inserted at the beginning, like the following
```
    BayesFactorBeta = 0.122298 *  w=0.124629.ctl
```
This specifies the beta value when the control file is used to run
BPP.

Second `BFdriver` creates a file named `betaweights.txt`, which lists
the beta values and Gauss-Legendre weights. I have copied those values
into an excel file in `frogs/BFdriver.frogs.xls`.

Third `BFdriver` creates a file named commands, which has the bash shell
scripts for submitting the 16 jobs using `qsub`. You use the following
to submit the jobs.
```
    source commands
```
You can look at the content of `tree1.sh` to see the script for the last
job:
```
    more tree1.sh
```
which should have the content like the following:
```
    bpp --cfile bpp.b16.ctl > log.b16.txt
```
You can use `qstat` to check the status of the 16 jobs you have
submitted. When the jobs are running, they generate output files in the
current folder, such as `mcmc.b01.txt`, `out.b01.txt`, and `log.b01.txt`
(which logs the screen output). After all jobs are finished, you can use
grep to extract the line with `BFbeta` from screen log (this command is
at the bottom of the file commands).
```
    grep BFbeta log.b*.txt
```
Then copy the ElnfX values into the excel file, and estimate the
logarithm of the marginal likelihood by summing (`weights * ElnfX / 2`)
over the 16 points. This gives the log marginal likelihood to be
$-3185.93$.

#### Exercise and results

Duplicate the calculation for the alternative species tree (tree2) in
the folder `/bf2`. Change the species tree topology to (((L, H), C), K),
and use tree2.sh as the temporary job script file. Everything else
should be the same as described above.

My runs gave the log marginal likelihood for tree 2 to be $-3186.14$.
The ratio of the posterior probabilities for the two trees is then
estimated to be

$\frac{P_1}{P_2} = \mathrm{e}^{-3185.93 - (-3186.14)} = \mathrm{e}^{0.21} = 1.24$. 

With BPP4.0 and the inverse-gamma priors on $\theta$ and
$\tau$s, the MCMC runs (A01) gave the posterior probabilities for trees
1 and 2 as 0.167 and 0.137, with the ratio 1.22. The MCMC runs and the
marginal likelihood calculations seem to agree with each other. (Note
that with BPP3 and the gamma priors on $\theta$ and
$\tau$s, the posteriors are 0.16 and 0.13 with the ratio 1.2 
([Yang 2015](https://doi.org/10.1093/czoolo/61.5.854), fig.4).

#### Common errors and problems

Check the bpp control file by running bpp at the command line before
submitting the jobs. Make sure that the bpp program is on your path or
use a full path. You may have to edit the source file `BFdriver.c` and
recompile.
