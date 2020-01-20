import lapper, tables, hts, strutils, re, algorithm

type
  Region* = ref object of RootObj
    chrom*: string
    start*: int
    stop*: int

  GFFGene* = ref object of Region
    gene_symbol: string
    gene_id: int
    biotype: string
    exons: seq[Region]

proc start*(region : Region): int {.inline.} = return region.start
proc stop*(region : Region): int {.inline.} = return region.stop
proc len*(region : Region): int {.inline.} = return region.stop - region.start + 1
proc `==`*(a, b: Region): bool {.inline.} = return a.chrom == b.chrom and a.start == b.start and a.stop == b.stop
proc merge*(a: var Region, b: Region) {.inline.} = a.start = min(a.start, b.start); a.stop = max(a.stop, b.stop)
proc globalFractionOverlap*(a: Region, b: Region): float {.inline.} = return (min(a.stop, b.stop) - max(a.start, b.start) + 1) / (max(a.stop, b.stop) - min(a.start, b.start) + 1)
proc fractionOverlap*(a: Region, b: Region): float {.inline.} = return (min(a.stop, b.stop) - max(a.start, b.start) + 1) / a.len()

# Coversion table from NCBI chromosomes ID's to usual names
let
  # FIXME: Should we remove the version from NC_ ids ?
  ncbi_to_chr = {
    "NC_000001.10": "1",
    "NC_000002.11": "2",
    "NC_000003.11": "3",
    "NC_000004.11": "4",
    "NC_000005.9": "5",
    "NC_000006.11": "6",
    "NC_000007.13": "7",
    "NC_000008.10": "8",
    "NC_000009.11": "9",
    "NC_000010.10": "10",
    "NC_000011.9": "11",
    "NC_000012.11": "12",
    "NC_000013.10": "13",
    "NC_000014.8": "14",
    "NC_000015.9": "15",
    "NC_000016.9": "16",
    "NC_000017.10": "17",
    "NC_000018.9": "18",
    "NC_000019.9": "19",
    "NC_000020.10": "20",
    "NC_000021.8": "21",
    "NC_000022.10": "22",
    "NC_000023.10": "X",
    "NC_000024.9": "Y",
    "NC_012920.1": "M"
    }.toTable

proc removeChrPrevix*(chrom: string): string =
  if chrom =~ re"""^chr(.*)""":
    return matches[0]
  else:
    return chrom

proc parseChr*(chrom: string): string {.inline.} =
  if ncbi_to_chr.hasKey(chrom):
    result = ncbi_to_chr[chrom]
  else:
    result = removeChrPrevix(chrom)

proc parseKeyValues*(str: string, global_sep: char, key_value_sep: char): TableRef[string, string] =
  let fields = str.split(global_sep)
  result = newTable[string, string]()
  for f in fields:
    let 
      kv_split = f.split(key_value_sep, 1)
    if kv_split.len() == 2:
      result[kv_split[0]] = kv_split[1]
    else:
      stderr.writeLine("[Error] Value fields " & f & " was not a key/value field using separator " & key_value_sep)

proc loadGenesFromGFF*(gff_file: string): TableRef[string, Lapper[GFFGene]] =
  result = newTable[string, Lapper[GFFGene]]() 
  var
    fh: BGZ
    genes_chr_table = newTable[string, seq[GFFGene]]() # Temp table to load genes per-chromosomes
    genes_name_table = newTable[string, GFFGene]() 

  open(fh, gff_file, "r")
  for line in fh:
    # Skip headers
    if line.len() == 0 or line[0] == '#':
      continue
    var v = line.split('\t', 3)
    
    # Only use "BestRefSeq" annotations
    if v[1] != "BestRefSeq":
      continue

    # NC_000001.10    BestRefSeq      gene    367659  368597  .       +       .       ID=gene-OR4F29;Dbxref=GeneID:729759,HGNC:HGNC:31275;Name=OR4F29;description=olfactory receptor family 4 subfamily F member 29;gbkey=Gene;gene=OR4F29;gene_biotype=protein_coding;gene_synonym=OR7-21
    if v[2] == "gene":
      var
        v2 = v[3].split('\t')
        gene = GFFGene(chrom: parseChr(v[0]), start: parseInt(v2[0]), stop: parseInt(v2[1]), exons: @[])
        gff_fields = v2[5].parseKeyValues(';','=')
        dbxref_fields = gff_fields["Dbxref"].parseKeyValues(',',':')
      
      gene.gene_symbol = gff_fields["Name"]
      if dbxref_fields.hasKey("GeneID"):
        gene.gene_id = parseInt(dbxref_fields["GeneID"])
      
      if gff_fields.hasKey("gene_biotype"):
        gene.biotype = gff_fields["gene_biotype"]
      
      if genes_chr_table.hasKey(gene.chrom):
        genes_chr_table[gene.chrom].add(gene)
      else:
        genes_chr_table[gene.chrom] = @[gene]
      
      genes_name_table[gene.gene_symbol] = gene

    # NC_000001.10    Curated Genomic exon    131068  132927  .       +       .       ID=id-CICP27;Parent=gene-CICP27;Dbxref=GeneID:100420257,HGNC:HGNC:48835;gbkey=exon;gene=CICP27
    elif v[2] == "exon":
      var
        v2 = v[3].split('\t')
        gff_fields = v2[5].parseKeyValues(';','=')
        dbxref_fields = gff_fields["Dbxref"].parseKeyValues(',',':')
        gene_symbol : string

      if gff_fields.hasKey("gene"):
        gene_symbol = gff_fields["gene"]
      
      if dbxref_fields.hasKey("Genbank"):
        if dbxref_fields["Genbank"].substr(0, 1) != "NM":
          continue

      # This exon belongs to an gene we are annotation, we catch it
      if gene_symbol != "" and genes_name_table.hasKey(gene_symbol):
        let
          exon = Region(chrom: parseChr(v[0]),start: parseInt(v2[0]), stop: parseInt(v2[1]))
        
        # Only add uniq exons and merge overlapping ones
        var i = 0
        for e in genes_name_table[gene_symbol].exons.mitems():
          let overlap_fraction = e.fraction_overlap(exon)
          if overlap_fraction == 1:
            break
          elif overlap_fraction > 0:
            e.merge(exon)
            break
          inc(i)
        
        # The exon has not been found / merge, we add it
        if i == genes_name_table[gene_symbol].exons.len():
          genes_name_table[gene_symbol].exons.add(exon)
      
  # TODO: Print number of loaded genes !!!
  
  # Load set of genes (per chromosome) to lapper index
  stderr.writeLine("[Log] Create lapper index for file " & gff_file)
  for chrom in genes_chr_table.keys():
    result[chrom] = lapify(genes_chr_table[chrom])

proc cmpGenes*(x, y: GFFGene): int =
  ## TODO: We should select Exonic over non-exonic regions (even if there is a non-coding gene)
  ## We select protein coding over non-coding gene (always ?)
  if x.biotype == "protein_coding" and y.biotype != "protein_coding":
    return -1
  elif x.biotype != "protein_coding" and y.biotype == "protein_coding":
    return 1
  else:
    # We select the oldest gene_id
    return cmp(x.gene_id, y.gene_id)

proc getInfoString*(genes_index: TableRef[string, Lapper[GFFGene]], chrom: string, start: int, stop: int): string =
  if genes_index.hasKey(chrom):
    var 
      res = new_seq[GFFGene]() # Store retrieved genes 
      found_overlapping_genes = genes_index[chrom].find(start, stop, res)

    if found_overlapping_genes:
      # Sort genes
      res.sort(cmpGenes)
      var gene_info: seq[string]
      for gene in res:
        gene_info.add(gene.gene_symbol & ":" & $gene.gene_id)
      result = gene_info.join("|")
  else:
    stderr.writeLine("[Error] Chrom " & chrom & " not found in GFF annotations")