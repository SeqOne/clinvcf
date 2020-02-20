![ClinVCF-logo](clinvcf.png)

# ClinVCF

ClinVCF **generates a VCF file from a Clinvar Full Release** (XML format). It was first developped because we observed missing variants in VCF files provided by NCBI. We later extended its capabilities to provived enhanced Clinvar VCF files by :

- **Improving Clinvar classification and aggregation method** by [deciphering "conflicting intepretation" records](#clinicalsignificance-correction-module) where almost all submissions goes into the same direction.
- **Implementing a more robust [gene annotation module](#gene-annotation)** based NCBI GFF files.

ClinVCF is **developped in NimLang, is highly efficient*** (~ 5 minutes to generate the VCF from the XML) and supports GRCh37 and GRCh38 genomes builds.

**Table of content**
- [ClinVCF](#clinvcf)
  - [Quick start](#quick-start)
  - [Usage](#usage)
    - [Output format](#output-format)
  - [Methodology](#methodology)
    - [ClinicalSignificance correction module](#clinicalsignificance-correction-module)
    - [Gene annotation](#gene-annotation)

## Quick start

You need to have [nimlang installed](https://nim-lang.org/install_unix.html) to compile and install clinVCF.

```
# Git clone and install
git clone https://gitlab.seq.one/workset/clinvcf.git && cd clinvcf && nimble install

# Download (latest) Clinvar XML release
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz

# Download GFF for gene annotation (GRCh37 or 38)
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz

# Generate clinvar VCF (GRCh37 by default)
clinvcf ClinVarFullRelease_00-latest.xml.gz | bgzip -c > clinvar_GRCh37.vcf.gz
```

## Usage

```
Usage: clinvcf [options] <clinvar.xml.gz>

Options:
  --genome <version>              Genome assembly to use [default: GRCh37]

Gene annotation:
  --gff <file>                    NCBI GFF to annotate variations with genes
  --coding-first                  Give priority to coding gene in annotation (even if intronic and exonic for another gene)
  --gene-padding <int>            Padding to annotation upstream/downstream genes (not applied for MT) [default: 5000]
```

### Output format

ClinVCF generates a VCF with almost identical format as the original NCBI VCF.

However, not all VCF fields are currently support by ClinVCF (see table bellow), and
additionnal fields are provided.

VCF Info field | Format | Description | Example
---------------|--------|-------------|---------
**ALLELEID** | *Integer* | the ClinVar Allele ID | `1234`
**CLNREVSTAT** | *String* | [ClinVar review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/) for the Variation ID | `no_assertion_criteria_provided` 
**CLNSIG** | String | [Clinical significance](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/) for this single variant | `Pathogenic/Likely_Pathogenic`
**OLD_CLNSIG**  | String | Orignial Clinical significance if variant reclassified by clinVCF correction module | `Conflicting_interpretations_of_pathogenicity`
**CLNRECSTAT** | Integer | [3-levels stars confidence](#clinicalsignificance-correction-module) of Variant Alert! automatic reclassfication. | `3`
**GENEINFO** | String | Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (`:`) and each pair is delimited by a vertical bar (`\|`) | `FTCD:10841\|FTCD-AS1:100861507`
**MC** | String | comma separated list of molecular consequence in the form of Sequence Ontology `ID\|molecular_consequence` | `SO:0001583\|missense_variant`
**RS** | String | dbSNP ID (i.e. rs number) | `80358507`

## Methodology

### ClinicalSignificance correction module

According to the 1.5 * IQR method, we remove outliers submissions and reclassify conflicting status variants according to ClinVar policies. We apply a 3-level star metrics according to our reclassification confidence. 4 or more submission is needed. We only reclassify variants from `conflicting` status to `likely pathogenic` or `pathogenic` status. 

- ⭐ **(1 star)** : default
- ⭐⭐ **(2 stars)** : reclassification remains even if we add a virtual VUS submission
- ⭐⭐⭐ **(3 stars)** : 2 stars requirements and at least 1 pathogenic classification

### Gene annotation

1. **We load all genes from the input GFF** and add them to the index with a padding (5000bp by default and 2bp for MT genes), to annotate upstream / downstream variants.
2. **For each variant we query the gene index** and retrieve all overlapping genes.
3. **Overlapped genes are later prioritize** in the `GENEINFO` field with two different procedures (depending of clinVCF parameter)
     - If `--coding-first` option is activated :
       - We take coding genes over all other genes (except for MT genome)
       - If we have an equality we take exonic (+/-20bp padding) over intronic/intergenic candidates
       - If none are exonic, we take the gene with closest exon
       - If both are exonic, we take the oldest gene ID in NCBI Entrez database
    - Default procedure :
      - We take coding gene over all other genes (except for MT genome) if the variant is exonic (+/- 20bp)
      - If we have an equality we take exonic (+/-20bp padding) over intronic/intergenic candidates
      - If none are exonic, we take the gene with closest exon
      - If both are exonic, we take the oldest gene ID in NCBI Entrez database