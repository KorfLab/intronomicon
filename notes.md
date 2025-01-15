Notes
=====

SRA, GEO, PubMed


## GEO ##

gpl - platform (e.g. Illumina NovoSeq 6000 w/ C. elegans)
	gse - a study (e.g. a paper or papers)
		gsm - samples --> SRX
			some controls
				which can contain multiple sequencing runs
			some experimentals
				which can contain multiple sequencing runs


acc = a valid GEO accession i.e., gplxxx, gsmxxx or gsexxx
targ = self, gsm, gpl, gse or all
view = brief, quick, data or full
form = text, html or xml


## GEO Example GSE237802 ##

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237802&targ=self&view=brief&form=text

"The 18S rRNA Methyltransferase DIMT-1 Regulates Lifespan in the Germline Later in Life"

PubMed: 38798397, 38946979

GSM:
GSM7653725 Germline ribosome-bound mRNA, control RNAi, biological replicate 1
GSM7653726 Germline ribosome-bound mRNA, control RNAi, biological replicate 2
GSM7653727 Germline ribosome-bound mRNA, control RNAi, biological replicate 3
GSM7653728 Germline ribosome-bound mRNA, control RNAi, biological replicate 4
GSM7653729 Germline ribosome-bound mRNA, dimt-1 RNAi, biological replicate 1
GSM7653730 Germline ribosome-bound mRNA, dimt-1 RNAi, biological replicate 2
GSM7653731 Germline ribosome-bound mRNA, dimt-1 RNAi, biological replicate 3
GSM7653732 Germline ribosome-bound mRNA, dimt-1 RNAi, biological replicate 4
GSM7653733 Total RNA, control RNAi, biological replicate 1
GSM7653734 Total RNA, control RNAi, biological replicate 2
GSM7653735 Total RNA, control RNAi, biological replicate 3
GSM7653736 Total RNA, control RNAi, biological replicate 4
GSM7653737 Total RNA, dimt-1 RNAi, biological replicate 1
GSM7653738 Total RNA, dimt-1 RNAi, biological replicate 2
GSM7653739 Total RNA, dimt-1 RNAi, biological replicate 3
GSM7653740 Total RNA, dimt-1 RNAi, biological replicate 4


https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7653725&targ=self&view=brief&form=text


Each GSM links directly to a SRX

BioProject: PRJNA996637 links directly to GSE237802








The text form is very easily parsed

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7653734&targ=self&view=brief&form=text

targ=gsm and targ=gse are very different

difference between GSM and GSE?

GSE is the parent description
GSM are like runs
	but could some runs be mutant and some be WT?
	I think so

Are some of these WT and some mutant?

!Series_sample_id = GSM7653725
!Series_sample_id = GSM7653726
!Series_sample_id = GSM7653727
!Series_sample_id = GSM7653728
!Series_sample_id = GSM7653729
!Series_sample_id = GSM7653730
!Series_sample_id = GSM7653731
!Series_sample_id = GSM7653732
!Series_sample_id = GSM7653733
!Series_sample_id = GSM7653734
!Series_sample_id = GSM7653735
!Series_sample_id = GSM7653736
!Series_sample_id = GSM7653737
!Series_sample_id = GSM7653738
!Series_sample_id = GSM7653739
!Series_sample_id = GSM7653740



!Series_title = The 18S rRNA Methyltransferase DIMT-1 Regulates Lifespan in the G
ermline Later in Life

!Series_pubmed_id = 38798397
!Series_pubmed_id = 38946979


!Series_summary = There is a marked dysregulation of the proteome as organisms ag
e. Regulation of both protein production and degradation can become dysregulated
as an organism ages. However, whether and how aging-responsive transcripts are se
lectively translated is unknown. Changes in the composition of the ribosome can l
ead to differential binding and translation of specific transcripts. Here we exam
ined the role of ribosomal RNA (rRNA) methylation in maintaining appropriate tran
slation of age-dependent transcripts. In a directed RNAi screen we identified the
 18S rRNA N6’-dimethyl adenosine (m6,2A) methyltransferase, dimt-1, as a regulato
r of C. elegans lifespan and stress resistance. Lifespan extension induced by dim
t-1 deficiency required a functional germline and the Rag GTPase, raga-1, which l
inks amino acid sensing to mechanistic target of rapamycin (mTOR) complex (mTORC)
1. Using an auxin-inducible degron tagged version of dimt-1 we demonstrate that D
IMT-1 functions in the worms germline after mid-life to regulate lifespan. We fur
ther found that knock-down of dimt-1 leads to selective translation of transcript
s important for stress resistance and lifespan regulation in the C. elegans germl
ine in mid-life. Our findings highlight a new role for rRNA modifications such as
 m6,2A and ribosome heterogeneity more broadly, in maintaining appropriate transl
ation later in life to promote healthy aging.

!Series_overall_design = To investigate the role of germline dimt-1 in regulating C. elegans lifespan, we IP'd FLAG-tagged germline ribosomes from dimt-1 knock down and control worms. The purified ribosome-bound transcripts were then subject to mRNA-seq. As a control, we also performed mRNA-seq on dimt-1 knock down and control worms using total RNA from whole-worm lysates.

!Series_platform_id = GPL26672
!Series_platform_organism = Caenorhabditis elegans
!Series_platform_taxid = 6239
!Series_sample_organism = Caenorhabditis elegans
!Series_sample_taxid = 6239
