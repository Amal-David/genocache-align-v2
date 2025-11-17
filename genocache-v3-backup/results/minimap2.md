Sequence analysis
Minimap2: pairwise alignment for nucleotide
sequences
Heng Li*
Department of Medical Population Genetics Program, Broad Institute, Cambridge, MA 02142, USA
*To whom correspondence should be addressed.
Associate Editor: Inanc Birol
Received on January 1, 2018; revised on March 16, 2018; editorial decision on March 22, 2018; accepted on May 4, 2018
Abstract
Motivation: Recent advances in sequencing technologies promise ultra-long reads of 100 kb in
average, full-length mRNA or cDNA reads in high throughput and genomic contigs over 100 Mb
in length. Existing alignment programs are unable or inefficient to process such data at scale,
which presses for the development of new alignment algorithms.
Results: Minimap2 is a general-purpose alignment program to map DNA or long mRNA
sequences against a large reference database. It works with accurate short reads of 100 bp in
length, 1 kb genomic reads at error rate 15%, full-length noisy Direct RNA or cDNA reads
and assembly contigs or closely related full chromosomes of hundreds of megabases in length.
Minimap2 does split-read alignment, employs concave gap cost for long insertions and dele-
tions and introduces new heuristics to reduce spurious alignments. It is 3–4 times as fast as
mainstream short-read mappers at comparable accuracy, and is 30 times faster than long-
read genomic or cDNA mappers at higher accuracy, surpassing most aligners specialized in
one type of alignment.
Availability and implementation: https://github.com/lh3/minimap2
Contact: hengli@broadinstitute.org
Supplementary information: Supplementary data are available at Bioinformatics online.
1 Introduction
Single Molecule Real-Time (SMRT) sequencing technology and
Oxford Nanopore technologies (ONT) produce reads over 10 kbp
in length at an error rate 15%. Several aligners have been devel-
oped for such data (Chaisson and Tesler, 2012; Li, 2013; Lin and
Hsu, 2017; Liu et al., 2016, 2017; Sedlazeck et al., 2018; Sovic
et al., 2016). Most of them were five times as slow as mainstream
short-read aligners (Langmead and Salzberg, 2012; Li, 2013) in
terms of the number of bases mapped per second. We speculated
there could be substantial room for speedup on the thought that
10 kb long sequences should be easier to map than 100 bp reads
because we can more effectively skip repetitive regions, which are
often the bottleneck of short-read alignment. We confirmed our
speculation by achieving approximate mapping 50 times faster than
BWA-MEM (Li, 2016). Suzuki and Kasahara (2018) extended our
work with a fast and novel algorithm on generating base-level
alignment, which in turn inspired us to develop minimap2 with
added functionality.
Both SMRT and ONT have been applied to the sequencing of
spliced mRNAs (RNA-seq). While traditional mRNA aligners work
(Iwata and Gotoh, 2012; Wu and Watanabe, 2005), they are not opti-
mized for long noisy sequence reads and are tens of times slower than
dedicated long-read aligners. When developing minimap2 initially for
aligning genomic DNA only, we realized minor modifications could
enable the base algorithm to map mRNAs as well. Minimap2
becomes a first RNA-seq aligner specifically designed for long noisy
reads. We have also extended the original algorithm to map short
reads at a speed faster than several mainstream short-read mappers.
In this article, we will describe the minimap2 algorithm and its
applications to different types of input sequences. We will evaluate
the performance and accuracy of minimap2 on several simulated
and real datasets and demonstrate the versatility of minimap2.
V C The Author(s) 2018. Published by Oxford University Press. All rights reserved. For permissions, please e-mail: journals.permissions@oup.com 3094
Downloaded from https://academic.oup.com/bioinformatics/article/34/18/3094/4994778 by guest on 13 November 2025
Aligning nucleotide sequences with minimap2 2 Materials and methods
Minimap2 follows a typical seed-chain-align procedure as is used
by most full-genome aligners. It collects minimizers (Roberts et al.,
2004) of the reference sequences and indexes them in a hash table,
with the key being the hash of a minimizer and the value being a
list of locations of the minimizer copies. Then for each query
sequence, minimap2 takes query minimizers as seeds, finds exact
matches (i.e. anchors) to the reference, and identifies sets of coli-
near anchors as chains. If base-level alignment is requested, mini-
map2 applies dynamic programming (DP) to extend from the ends
of chains and to close regions between adjacent anchors in chains.
Minimap2 uses indexing and seeding algorithms similar to mini-
map (Li, 2016), and furthers the predecessor with more accurate
chaining, the ability to produce base-level alignment and the support
of spliced alignment.
2.1 Chaining
2.1.1 Finding the optimal chaining scores
An anchor is a 3-tuple (x, y, w), indicating interval x w þ1; x ½ on
the reference matching interval y w þ1; y ½ on the query. Given a
list of anchors sorted by ending reference position x, let f(i) be the
maximal chaining score up to the i-th anchor in the list. f(i) can be
calculated with dynamic programming:
f i ðÞ¼maxfmax
i>j 1ff j ðÞþa j; i ð Þb j; i ð Þg; wig (1)
where a j; i ð Þ¼minfminfyi yj; xi xjg; wig is the number of
matching bases between the two anchors. b j; i ð Þ> 0 is the gap cost.
It equals 1if yj yi or maxfyi yj; xi xjg> G (i.e. the distance
between two anchors is too large); otherwise
b j; i ð Þ¼cc yi yj xi xj (2)
In implementation, a gap of length l costs
cc l ðÞ¼
(
0:01 w jljþ0:5 log2jlj l¼0 ð Þ
0 l¼0 ð Þ
where w is the average seed length. For N anchors, directly comput-
ing all f ðÞwith Equation (1) takes O N2 time. Although theoreti-
cally faster chaining algorithms exist (Abouelhoda and Ohlebusch,
2005), they are inapplicable to generic gap cost, complex to imple-
ment and usually associated with a large constant. We introduced a
simple heuristic to accelerate chaining.
We note that if anchor i is chained to j, chaining i to a predecessor
of j is likely to yield a lower score. When evaluating Equation (1), we
start from anchor i 1 and stop the process if we cannot find a better
score after up to h iterations. This approach reduces the average time
to O(hN). In practice, we can almost always find the optimal chain
with h¼50; even if the heuristic fails, the optimal chain is often close.
2.1.2 Backtracking
Let P(i) be the index of the best predecessor of anchor i. It equals 0
if f i ðÞ¼wi or argmaxjff j ðÞþa j; i ð Þb j; i ð Þgotherwise. For each
anchor i in the descending order of f(i), we apply PðÞrepeatedly to
find its predecessor and mark each visited i as ‘used’, until P(i)¼0
or we reach an already ‘used’ i. This way we find all chains with no
anchors used in more than one chains.
2.1.3 Identifying primary chains
In the absence of copy number changes, each query segment should
not be mapped to two places in the reference. However, chains
3095
found at the previous step may have significant or complete overlaps
due to repeats in the reference (Li and Durbin, 2010). Minimap2
used the following procedure to identify primary chains that do not
greatly overlap on the query.
Let Q be an empty set initially. For each chain from the best to the
worst according to their chaining scores: if on the query, the chain
overlaps with a chain in Q by 50% or higher percentage of the shorter
chain, mark the chain as secondary to the chain in Q; otherwise, add
the chain to Q. In the end, Q contains all the primary chains. We did
not choose a more sophisticated data structure (e.g. range tree or k-d
tree) because this step is not the performance bottleneck.
For each primary chain, minimap2 estimates its mapping quality
with an empirical formula:
mapQ¼40 1 f2=f1 ð Þminf1; m=10glog f1
where log denotes natural logarithm, m is the number of anchors on
the primary chain, f1 is the chaining score, and f2 f1 is the score
of the best chain that is secondary to the primary chain. Intuitively,
a chain is assigned to a higher mapping quality if it is long and its
best secondary chain is weak.
2.1.4 Estimating per-base sequence divergence
Suppose a query sequence harbors n seeds of length k, m of which
are present in a chain. We want to estimate the sequence divergence
between the query and the reference sequences in the chain. This is
useful when base-level alignment is too expensive to perform.
If we model substitutions with a homogeneous Poisson process
along the query sequence, the probablity of seeing k consecutive
bases without substitutions is e k . On the assumption that all
k-mers are independent of each other, the likelihood function of is
L jn; m; k ð Þ¼e m k 1 e k n m
The maximum likelihood estimate of is
b
¼
1
k log n
m:
In reality, sequencing errors are sometimes clustered and k-mers are
not independent of each other, especially when we take minimizers as
seeds. These violate the assumptions in the derivation above. As a
result,b is only approximate and can be biased. It also ignores long
deletions from the reference sequence. In practice, fortunately,b is
often close to and strongly correlated with the sequence divergence
estimated from base-level alignments. On the several datasets used in
Section 3.1, the Spearman correlation coefficient is around 0.9.
2.1.5 Indexing with homopolymer compressed k-mers
SmartDenovo (https://github.com/ruanjue/smartdenovo; J. Ruan,
personal communication) indexes reads with homopolymer-
compressed (HPC) k-mers and finds the strategy improves overlap
sensitivity for SMRT reads. Minimap2 adopts the same heuristic.
The HPC string of a string s, denoted by HPC s ðÞ, is constructed by
contracting homopolymers in s to a single base. An HPC k-mer of s is a
k-long substring of HPC s ðÞ. For example, suppose s ¼GGATTTTCCA,
HPC s ðÞ¼GATCA and the first HPC 4-mer is GATC.
To demonstrate the effectiveness of HPC k-mers, we performed
read overlapping for the example E. coli SMRT reads from PBcR
(Berlin et al., 2015), using different types of k-mers. With normal
15 bp minimizers per 5 bp window, minimap2 finds 90.9% of 2 kb
overlaps inferred from the read-to-reference alignment. With HPC
19-mers per 5 bp window, minimap2 finds 97.4% of overlaps.
Downloaded from https://academic.oup.com/bioinformatics/article/34/18/3094/4994778 by guest on 13 November 2025
3096 H.Li
It achieves this higher sensitivity by indexing 1/3 fewer minimizers,
which further helps performance. HPC-based indexing reduces the
sensitivity for current ONT reads, though.
2.2 Aligning genomic DNA
2.2.1 Alignment with 2-piece affine gap cost
Minimap2 performs DP-based global alignment between adjacent
anchors in a chain. It uses a 2-piece affine gap cost (Gotoh, 1990):
ca l ðÞ¼minfq þjlje;~
q þjlj~
eg (3)
Without losing generality, we always assume q þe <~
q þ~
e. On the
condition that e >~
e, it applies cost q þjlje to gaps shorter than d
~
q q ð Þ= e~
e ð Þeand applies~
q þjlj~
e to longer gaps. This scheme
helps to recover longer insertions and deletions (INDELs).
The equation to compute the optimal alignment under caðÞis
8 > > > > > > > > > < > > > > > > > > > :
Hij¼maxfHi 1;j 1 þs i; j ð Þ; Eij; Fij ;
Eiþ1;j¼maxfHij q; Eijge
Fi;jþ1 ¼maxfHij q; Fijge
~
~
Eiþ1;j¼maxfHij
~
q;
e
Eijg~
~
~
Fi;jþ1 ¼maxfHij
~
q;
e
Fijg~
~
Eij;
~
Fijg
(4)
where s(i, j) is the score between the i-th reference base and j-th
query base. Equation (4) is a natural extension to the equation under
affine gap cost (Altschul and Erickson, 1986; Gotoh, 1982).
2.2.2 The Suzuki–Kasahara formulation
When we allow gaps longer than several hundred base pairs,
nucleotide-level alignment is much slower than chaining. SSE accel-
eration is critical to the performance of minimap2. Traditional SSE
implementations (Farrar, 2007) based on Equation (4) can achieve
16-way parallelization for short sequences, but only 4-way paralleli-
zation when the peak alignment score reaches 32 767. Long
sequence alignment may exceed this threshold. Inspired by Wu et al.
(1996) and the following work, Suzuki and Kasahara (2018) pro-
posed a difference-based formulation that lifted this limitation. In
case of 2-piece gap cost, define
8 > > < > > :
uij¢Hij Hi 1;j vij¢Hij Hi;j 1
xij¢Eiþ1;j Hij
~
xij¢~
Eiþ1;j Hij
yij¢Fi;jþ1 Hij
~
yij¢~
Fi;jþ1 Hij
We can transform Equation (4) to
8 > > > > > > > > > > > > > > > < > > > > > > > > > > > > > > > > > > > :
zij¼ maxfs i; j ð Þ; xi 1;j þvi 1;j; yi;j 1 þui;j 1;
~
xi 1;j þvi 1;j;~
yi;j 1 þui;j 1g
uij¼ zij vi 1;j
vij¼ zij ui;j 1
xij¼ maxf0; xi 1;j þvi 1;j zij þqgq e
yij¼ maxf0; yi;j 1 þui;j 1 zij þqgq e
~
xij¼ maxf0;~
xi 1;j þvi 1;j zij þ~
qg~
q~
e
~
yij¼ maxf0;~
yi;j 1 þui;j 1 zij þ~
qg~
q~
e
(5)
where zij is a temporary variable that does not need to be stored.
An important property of Equation (5) is that all values are
bounded by scoring parameters. To see that,
xij¼Eiþ1;j Hij¼maxfq; Eij Hij ge
With Eij Hij, we have
q e xij maxfq; 0ge ¼e
and similar inequations for yij,~
xij and~
yij. In addition,
uij¼zij vi 1;j maxfxi 1;j;~
xi 1;jgq e
As the maximum value of zij¼Hij Hi 1;j 1 is M, the maximal
matching score, we can derive
uij M vi 1;j M þq þe
In conclusion, in Equation (5), x and y are bounded by q e; e ½ ,
~
x and~
y by~
q ½ ~
e;~
e , and u and v by q e; M þq þe ½ . When
128 q e < M þq þe 127, each of them can be stored
as a 8-bit integer. This enables 16-way SSE vectorization regardless
of the peak score of the alignment.
For a more efficient SSE implementation, we transform the row-
column coordinate to the diagonal-antidiagonal coordinate by let-
ting r i þj and t i. Equation (5) becomes:
8 > > > > > > > > > > > > > > > > > < > > > > > > > > > > > > > > > > > :
zrt ¼ maxfs t; r t ð Þ; xr 1;t 1 þvr 1;t 1; yr 1;t
þur 1;t ;~
xr 1;t 1 þvr 1;t 1;~
yr 1;t þur 1;t g
urt ¼ zrt vr 1;t 1
vrt ¼ zrt ur 1;t
xrt ¼ maxf0; xr 1;t 1 þvr 1;t 1 zrt þqgq e
yrt ¼ maxf0; yr 1;t þur 1;t zrt þqgq e
~
xrt ¼ maxf0;~
xr 1;t 1 þvr 1;t 1 zrt þ~
qg~
q~
e
~
yrt ¼ maxf0;~
yr 1;t þur 1;t zrt þ~
qg~
q~
e
In this formulation, cells with the same diagonal index r are inde-
pendent of each other. This allows us to fully vectorize the computa-
tion of all cells on the same anti-diagonal in one inner loop. It also
simplifies banded alignment (500 bp band width by default), which
would be difficult with striped vectorization (Farrar, 2007).
On the condition that q þe <~
q þ~
e and e >~
e, the initial values
in the diagonal-antidiagonal formuation are
8 > > < > > :
xr 1; 1 ¼yr 1;r ¼q e
~
xr 1; 1 ¼~
yr 1;r
¼~
q~
e
ur 1;r ¼vr 1; 1 ¼g r ðÞ
where
g r ðÞ¼
8 > > > > > > > > > > > < > > > > > > > > > r e~
q q ð Þ~
e r ¼~
q q
e ð Þ~
> > :
q e r ¼0 ð Þ
e r <
~
q q
1
e~
e
1
e~
e
~
e r >
~
q q
1
e~
e
These can be derived from the initial values for Equation (4).
When performing global alignment, we do not need to compute
Hrt in each cell. We use 16-way vectorization throughout the align-
ment process. When extending alignments from ends of chains, we
need to find the cell (r, t) where Hrt reaches the maximum. We resort
to 4-way vectorization to compute Hrt ¼Hr 1;t þurt . Because this
computation is simple, Equation (5) is still the dominant perform-
ance bottleneck.
In practice, our 16-way vectorized implementation of global
alignment is three times as fast as Parasail’s 4-way vectorization
Downloaded from https://academic.oup.com/bioinformatics/article/34/18/3094/4994778 by guest on 13 November 2025
Aligning nucleotide sequences with minimap2 (Daily, 2016). Without banding, our implementation is slower than
Edlib (Sosic and Sikic, 2017), but with a 1000 bp band, it is consid-
erably faster. When performing global alignment between anchors,
we expect the alignment to stay close to the diagonal of the DP
matrix. Banding is applicable most of the time.
2.2.3 The Z-drop heuristic
With global alignment, minimap2 may force to align unrelated
sequences between two adjacent anchors. To avoid such an artifact,
we compute accumulative alignment score along the alignment path
and break the alignment where the score drops too fast in the diago-
nal direction. More precisely, let S(i, j) be the alignment score along
the alignment path ending at cell (i, j) in the DP matrix. We break the
alignment if there exist i0
; j0 ð Þand (i, j), i0 < i and j0 < j, such that
S i0
; j0 ð ÞS i; j ð Þ> Z þe ji i0 ð Þ j j0 ð Þj
where e is the gap extension cost and Z is an arbitrary threshold.
This strategy is first used in BWA-MEM. It is similar to X-drop
employed in BLAST (Altschul et al., 1997), but unlike X-drop, it
would not break the alignment in the presence of a single long gap.
When minimap2 breaks a global alignment between two
anchors, it performs local alignment between the two subsequences
involved in the global alignment, but this time with the one subse-
quence reverse complemented. This additional alignment step may
identify short inversions that are missed during chaining.
2.2.4 Filtering out misplaced anchors
Due to sequencing errors and local homology, some anchors in a
chain may be wrong. If we blindly align regions between two mis-
placed anchors, we will produce a suboptimal alignment. To reduce
this artifact, we filter out anchors that lead to a >10 bp insertion
and a >10 bp deletion at the same time, and filter out terminal
anchors that lead to a long gap towards the ends of a chain. These
heuristics greatly alleviate the issues with misplaced anchors, but
they are unable to fix all such errors. Local misalignment is a limita-
tion of minimap2 which we hope to address in future.
2.3 Aligning spliced sequences
The algorithm described above can be adapted to spliced alignment.
In this mode, the chaining gap cost distinguishes insertions to and
deletions from the reference: cc l ðÞin Equation (2) takes the form of
cc l ðÞ¼
(
0:01 w l þ0:5 log2l l > 0 ð Þ
minf0:01 w jlj; log2jljg l < 0 ð Þ
Similarly, the gap cost function used for DP-based alignment is
changed to
ca l ðÞ¼ q þl e l > 0 ð Þ
(
minfq þjlje;~
qg l < 0 ð Þ
In alignment, a deletion no shorter than d~
q q ð Þ=eeis regarded
as an intron, which pays no cost to gap extensions.
To pinpoint precise splicing junctions, minimap2 introduces
reference-dependent cost to penalize non-canonical splicing:
8 > > > > > < > > > > > :
Hij¼maxfHi 1;j 1 þs i; j ð Þ; Eij ; Fij;
Eiþ1;j¼maxfHij q; Eijge
Fi;jþ1 ¼maxfHij q; Fijge
~
Eiþ1;j¼maxfHij d i ðÞ~
q;
~
Eijg
~
Eij a i ðÞg
(6)
3097
Let T be the reference sequence. d(i) is computed as
d i ðÞ¼
8 > > < > > :
0 if T i þ1; i þ3 ½ is GTA or GTG
p=2 if T i þ1; i þ3 ½ is GTC or GTT
p otherwise
where T i; j ½ extracts a substring of T between i and j inclusively. d(i)
penalizes non-canonical donor sites with p and less frequent
Eukaryotic splicing signal GT C=T ½ with p=2 (Irimia and Roy,
2008). Similarly,
a i ðÞ¼
8 > > < > > :
0 if T i 2; i ½ is CAG or TAG
p=2 if T i 2; i ½ is AAG or GAG
p otherwise
models the acceptor signal. Equation (6) is close to an equation in
Zhang and Gish (2006) except that we allow insertions immediately
followed by deletions and vice versa; in addition, we use the Suzuki–
Kasahara diagonal formulation in actual implementation.
If RNA-seq reads are not sequenced from stranded libraries, the
read strand relative to the underlying transcript is unknown. By
default, minimap2 aligns each chain twice, first assuming GT–AG as
the splicing signal and then assuming CT–AC, the reverse comple-
ment of GT–AG, as the splicing signal. The alignment with a higher
score is taken as the final alignment. This procedure also infers the
relative strand of reads that span canonical splicing sites.
In the spliced alignment mode, minimap2 further increases the
density of minimizers and disables banded alignment. Together with
the two-round DP-based alignment, spliced alignment is several
times slower than genomic DNA alignment.
2.4 Aligning short paired-end reads
During chaining, minimap2 takes a pair of reads as one fragment
with a gap of unknown length in the middle. It applies a normal gap
cost between seeds on the same read but is a more permissive gap
cost between seeds on different reads. More precisely, the gap cost
during chaining is (l¼0):
cc l ðÞ¼
(
0:01 w jljþ0:5 log2jlj if two seeds on the same read
minf0:01 w jlj; log2jljg otherwise
After identifying primary chains (Section 2.1.3), we split each frag-
ment chain into two read chains and perform alignment for each
read as in Section 2.2. Finally, we pair hits of each read end to find
consistent paired-end alignments.
3 Results
Minimap2 is implemented in the C programming language and
comes with APIs in both C and Python. It is distributed under the
MIT license, free to both commercial and academic uses. Minimap2
uses the same base algorithm for all applications, but it has to apply
different sets of parameters depending on input data types. Similar
to BWA-MEM, minimap2 introduces ‘presets’ that modify multiple
parameters with a simple invocation. Detailed settings and
command-line options can be found in the minimap2 manpage (sup-
plementary information). In addition to the applications evaluated
in the following sections, minimap2 also retains minimap’s function-
ality to find overlaps between long reads and to search against large
multi-species databases such as nt from NCBI.
Downloaded from https://academic.oup.com/bioinformatics/article/34/18/3094/4994778 by guest on 13 November 2025
3098 H.Li
3.1 Aligning long genomic reads
As a sanity check, we evaluated minimap2 on simulated human
reads along with BLASR (v1.MC.rc64; Chaisson and Tesler, 2012),
BWA-MEM (v0.7.15; Li, 2013), GraphMap (v0.5.2; Sovic et al.,
2016), Kart (v2.2.5; Lin and Hsu, 2017), minialign (v0.5.3; https://
github.com/ocxtal/minialign) and NGMLR (v0.2.5; Sedlazeck et al.,
2018). We excluded rHAT (Liu et al., 2016) and LAMSA (Liu et al.,
2017) because they either crashed or produced malformatted out-
put. In this evaluation, minimap2 has higher power to distinguish
unique and repetitive hits, and achieves overall higher mapping
accuracy (Fig. 1a). Minimap2 and NGMLR provide better mapping
quality estimate: they rarely give repetitive hits high mapping qual-
ity. Apparently, other aligners may occasionally miss close subopti-
mal hits and be overconfident in wrong mappings. On run time,
minimap2 took 200 CPU seconds, comparable to minialign and
Kart, and is over 30 times faster than the rest. Minimap2 consumed
6.8GB memory at the peak, more than BWA-MEM (5.4GB), similar
to NGMLR and less than others.
On real human SMRT reads, the relative performance and frac-
tion of mapped reads reported by these aligners are broadly similar
to the metrics on simulated data. We are unable to provide a good
estimate of mapping error rate due to the lack of the truth. On ONT
100 kb human reads (Jain et al., 2018), BWA-MEM failed. Kart,
minialign and minimap2 are over 70 times faster than others. We
have also examined tens of 100 bp INDELs in IGV (Robinson
et al., 2011) and can confirm the observation by Sedlazeck et al.
(2018) that BWA-MEM often breaks them into shorter gaps. The
issue is much alleviated with minimap2, thanks to the 2-piece affine
gap cost.
3.2 Aligning long spliced reads
We evaluated minimap2 on SIRV control data (AC: SRR5286959;
Byrne et al., 2017) where the truth is known. Minimap2 predicted
59 918 introns from 11 018 reads. 93.8% of splice juctions are pre-
cise. We examined wrongly predicted junctions and found the
majority were caused by clustered splicing signals (e.g. two adjacent
GT sites). When INDEL sequencing errors are frequent, it is difficult
to find precise splicing sites in this case. If we allow up to 10 bp dis-
tance from true splicing sites, 98.4% of aligned introns are approxi-
mately correct. It is worth noting that for SIRV, we asked minimap2
to model the GT::AG splicing signal only without extra bases. This
is because SIRV does not honor the evolutionarily prevalent signal
GT A=G ½ :: C=T ½ AG (Irimia and Roy, 2008).
We next aligned real mouse reads (Byrne et al., 2017) with
GMAP (v2017-06-20; Wu and Watanabe, 2005), minimap2, SpAln
(v2.3.1; Iwata and Gotoh, 2012) and STAR (v2.5.3a; Dobin et al.,
2013). In general, minimap2 is more consistent with existing anno-
tations (Table 1): it finds more junctions with a higher percentage
being exactly or approximately correct. Minimap2 is over 40 times
faster than GMAP and SpAln. While STAR is close to minimap2 in
speed, it does not work well with noisy reads.
We have also evaluated spliced aligners on a human Nanopore
Direct RNA-seq dataset (http://bit.ly/na12878ont). Minimap2
aligned 10 million reads in <1 wall-clock hour using 16 CPU cores.
94.2% of aligned splice junctions consistent with gene annotations.
In comparison, GMAP under option ‘-k 14 -n 0 -min-intronlength
30 –cross-species’ is 160 times slower; 68.7% of GMAP junctions
are found in known gene annotations. The percentage increases to
84.1% if an aligned junction within 10 bp from an annotated junc-
tion is considered to be correct. On a public Iso-Seq dataset (human
Alzheimer brain from http://bit.ly/isoseqpub), minimap2 is also
(a)
Fraction of mapped reads
(b)
Fraction of mapped reads
bowtie2
bwa-mem
minimap2
snap
1
1
0.98
0.98
0.96
0.96
0.94
0.92
0.9
blasr-mc
0.94
bwa-mem
graphmap
minialign
0.92
minimap2
ngm-lr 0.9
10-5 10-4 10-3 10-2 10-1
Error rate of mapped PacBio reads
10-7 10-6 10-5 10-4 10-3 10-2 10-1
Error rate of mapped short reads
Fig. 1. Evaluation on aligning simulated reads. Simulated reads were mapped
to the primary assembly of human genome GRCh38. A read is considered
correctly mapped if its longest alignment overlaps with the true interval, and
the overlap length is 10% of the true interval length. Read alignments are
sorted by mapping quality in the descending order. For each mapping quality
threshold, the fraction of alignments (out of the number of input reads) with
mapping quality above the threshold and their error rate are plotted along the
curve. (a) Long-read alignment evaluation. 33 088 1000 bp reads were simu-
lated using pbsim (Ono et al., 2013) with error profile sampled from file
‘m131017_060208_42213_*.1.*’ downloaded at http://bit.ly/chm1p5c3. The
N50 read length is 11 628. Aligners were run under the default setting for
SMRT reads. Kart outputted all alignments at mapping quality 60, so is not
shown in the figure. It mapped nearly all reads with 4.1% of alignments being
wrong, less accurate than others. (b) Short-read alignment evaluation. 10 mil-
lion pairs of 150 bp reads were simulated using mason2 (Holtgrewe, 2010)
with option ‘–illumina-prob-mismatch-scale 2.5’. Short-read aligners were
run under the default setting except for changing the maximum fragment
length to 800 bp
Table 1. Evaluation of junction accuracy on 2D ONT reads
GMAP minimap2 SpAln STAR
Run time (CPU min) 631 15.9 2076 33.9
Peak RAM (GByte) 8.9 14.5 3.2 29.2
# aligned reads 103 669 104 199 103 711 26 479
# chimeric alignments 1904 1488 0 0
# non-spliced alignments 15 854 14 798 17 033 10 545
# aligned introns 692 275 693 553 692 945 78 603
# novel introns 11 239 3113 8550 1214
% exact introns 83.8 94.0 87.9 55.2
% approx. introns 91.8 96.9 92.5 82.4
Notes: Mouse cDNA reads (AC: SRR5286960; R9.4 chemistry) were
mapped to the primary assembly of mouse genome GRCm38 with the follow-
ing tools and command options: minimap2 (‘-ax splice’); GMAP (‘-n 0 –min-
intronlength 30 –cross-species’); SpAln (‘-Q7 -LS -S3’); STARlong (according
to http://bit.ly/star-pb). The alignments were compared to the EnsEMBL gene
annotation, release 89. A predicted intron is novel if it has no overlaps with
any annotated introns. An intron is exact if it is identical to an annotated
intron. An intron is approximate if both its 50- and 30-end are within 10 bp
around the ends of an annotated intron. Chimeric alignments are defined in
the SAM spec (Li et al., 2009).
Downloaded from https://academic.oup.com/bioinformatics/article/34/18/3094/4994778 by guest on 13 November 2025
faster at higher junction accuracy in comparison to other aligners in
Table 1.
We noted that GMAP and SpAln have not been optimized for
noisy reads. We are showing the best setting we have experimented,
but their developers should be able to improve their accuracy
further.
3.3 Aligning short genomic reads
We evaluated minimap2 along with Bowtie2 [v2.3.3; (Langmead
and Salzberg, 2012)], BWA-MEM and SNAP [v1.0beta23; (Zaharia
Aligning nucleotide sequences with minimap2 et al., 2011)]. Minimap2 is 3–4 times as fast as Bowtie2 and BWA-
MEM, but is 1.3 times slower than SNAP. Minimap2 is more
accurate on this simulated dataset than Bowtie2 and SNAP but less
accurate than BWA-MEM (Fig. 1b). Closer investigation reveals
that BWA-MEM achieves a higher accuracy partly because it tries to
locally align a read in a small region close to its mate. If we disable
this feature, BWA-MEM becomes slightly less accurate than mini-
map2. We might implement a similar heuristic in minimap2 in
future.
To evaluate the accuracy of minimap2 on real data, we aligned
human reads (AC: ERR1341796) with BWA-MEM and minimap2,
and called SNPs and small INDELs with GATK HaplotypeCaller
v3.5 (Depristo et al., 2011). This run was sequenced from experi-
mentally mixed CHM1 and CHM13 cell lines. Both of them are
homozygous across the whole genome and have been de novo
assembled with SMRT reads to high quality. This allowed us to con-
struct an independent truth variant dataset (Li et al., 2017) for
ERR1341796. In this evaluation, minimap2 has higher SNP false
negative rate (FNR; 2.6% of minimap2 vs 2.3% of BWA-MEM),
but fewer false positive SNPs per million bases (FPPM; 7.0 versus
8.8), similar INDEL FNR (11.2% versus 11.3%) and similar
INDEL FPPM (6.4 versus 6.5). Minimap2 is broadly comparable to
BWA-MEM in the context of small variant calling.
3.4 Aligning long-read assemblies
Minimap2 can align a SMRT assembly (AC: GCA_001297185.1)
against GRCh38 in 7 min using eight CPU cores, over 20 times faster
than nucmer from MUMmer4 (Marc ¸ais et al., 2018). With the paf-
tools.js script from the minimap2 package, we called 2.67 million
single-base substitutions out of 2.78Gbp genomic regions. The
transition-to-transversion ratio (ts/tv) is 2.01. In comparison, using
MUMmer4’s dnadiff pipeline, we called 2.86 million substitutions
in 2.83Gbp at ts/tv¼1.87. Given that ts/tv averaged across the
human genome is about 2 but ts/tv averaged over random errors is
0.5, the minimap2 callset arguably has higher precision at lower
sensitivity.
The sample being assembled is a female. Minimap2 still called
201 substitutions on the Y chromosome. These substitutions all
come from one contig aligned at 96.8% sequence identity. The con-
tig could be a segmental duplication absent from GRCh38. In con-
strast, dnadiff called 9070 substitutions on the Y chromosome
across 73 SMRT contigs. This again implies our minimap2-based
pipeline has higher precision.
4 Discussions
Minimap2 is a versatile mapper and pairwise aligner for nucleotide
sequences. It works with short reads, assembly contigs and long
noisy genomic and RNA-seq reads, and can be used as a read map-
per, long-read overlapper or a full-genome aligner. Minimap2 is
also accurate and efficient, often outperforming other domain-
specific alignment tools in terms of both speed and accuracy.
The capability of minimap2 comes from a fast base-level align-
ment algorithm and an accurate chaining algorithm. When aligning
long query sequences, base-level alignment is often the performance
bottleneck. The Suzuki–Kasahara algorithm greatly alleviates the
bottleneck and enables DP-based splice alignment invol-
ving >100 kb introns, which was impractically slow 10 years ago.
The minimap2 chaining algorithm is fast and highly accurate by
itself. In fact, chaining alone is more accurate than all the other
long-read mappers in Figure 1a (data not shown). This accuracy
3099
helps to reduce downstream base-level alignment of candidate
chains, which is still several times slower than chaining even with
the Suzuki–Kasahara improvement. In addition, taking a general
form, minimap2 chaining can be adapted to non-typical data types
such as spliced reads and multiple reads per fragment. This gives us
the opportunity to extend the same base algorithm to a variety of
use cases.
Modern mainstream aligners often use a full-text index, such as
suffix array or FM-index, to index reference sequences. An advant-
age of this approach is that we can use exact seeds of arbitrary
lengths, which helps to increase seed uniqueness and reduce unsuc-
cessful extensions. Minimap2 indexes reference k-mers with a hash
table instead. Such fixed-length seeds are inferior to variable-length
seeds in theory, but can be computed much more efficiently in prac-
tice. When a query sequence has multiple seed hits, we can afford to
skip highly repetitive seeds without affecting the final accuracy. This
further alleviates the concern with the seeding uniqueness. At the
same time, at low sequence identity, it is rare to see long seeds any-
way. Hash table is the ideal data structure for mapping long noisy
sequences.
Acknowledgements
We owe a debt of gratitude to H. Suzuki and M. Kasahara for releasing their
masterpiece and insightful notes before formal publication. We thank M.
Schatz, P. Rescheneder and F. Sedlazeck for pointing out the limitation of
BWA-MEM. We are also grateful to minimap2 users who have greatly helped
to suggest features and to fix various issues.
Funding
This work was supported by the National Human Genome Research Institute
[1R01HG010040-01].
References
Abouelhoda,M.I. and Ohlebusch,E. (2005) Chaining algorithms for multiple
genome comparison. J. Discrete Algorithms, 3, 321–341.
Altschul,S.F. and Erickson,B.W. (1986) Optimal sequence alignment using
affine gap costs. Bull. Math. Biol., 48, 603–616.
Altschul,S.F. et al. (1997) Gapped BLAST and PSI-BLAST: a new generation
of protein database search programs. Nucleic Acids Res., 25, 3389–3402.
Berlin,K. et al. (2015) Assembling large genomes with single-molecule
sequencing and locality-sensitive hashing. Nat. Biotechnol., 33, 623–630.
Byrne,A. et al. (2017) Nanopore long-read RNAseq reveals widespread tran-
scriptional variation among the surface receptors of individual B cells. Nat.
Commun., 8, 16027.
Chaisson,M.J. and Tesler,G. (2012) Mapping single molecule sequencing
reads using basic local alignment with successive refinement (BLASR): appli-
cation and theory. BMC Bioinformatics, 13, 238.
Daily,J. (2016) Parasail: sIMD C library for global, semi-global, and local
pairwise sequence alignments. BMC Bioinformatics, 17, 81.
Depristo,M.A. et al. (2011) A framework for variation discovery and genotyp-
ing using next-generation DNA sequencing data. Nat. Genet., 43, 491–498.
Dobin,A. et al. (2013) STAR: ultrafast universal RNA-seq aligner.
Bioinformatics, 29, 15–21.
Farrar,M. (2007) Striped Smith-Waterman speeds database searches six times
over other SIMD implementations. Bioinformatics, 23, 156–161.
Gotoh,O. (1982) An improved algorithm for matching biological sequences.
J. Mol. Biol., 162, 705–708.
Gotoh,O. (1990) Optimal sequence alignment allowing for long gaps. Bull.
Math. Biol., 52, 359–373.
Holtgrewe,M. (2010) Mason: a read simulator for second generation sequenc-
ing data. Technical Report TR-B-10-06. Institut fu ¨r Mathematik und
Informatik, Freie Universita ¨t, Berlin