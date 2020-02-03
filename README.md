![ClinVCF-logo](clinvcf.png)

# ClinVCF or Clean Clinvar VCF

ClinVCF is a tool designed to generate a VCF from a Clinvar Full Release (XML format).

This tool has been developped as the VCF files provided by NCBI are missing variants, and
even though the error was reported, no improvement have been observed.

This tool is also designed to improve Clinvar classification and aggregation method as
deciphering "conflicting intepretation" for records where almost all submissions goes into the same direction.

This tool is also a way for SeqOne to provide an enhanced Clinvar VCF files that will be used in
SeqOne resources such has, i/ variant annotation, ii/ acmg classification module, iii/ variant alert notification.

## Usage

```
# Git clone and install
git clone https://gitlab.seq.one/workset/clinvcf.git && cd clinvcf && nimble install

# Download (latest) Clinvar XML release
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz

# Generate clinvar VCF (GRCh37 by default)
clinvcf ClinVarFullRelease_00-latest.xml.gz | bgzip -c > clinvar_GRCh37.vcf.gz
```

## ClinicalSignificance correction module

