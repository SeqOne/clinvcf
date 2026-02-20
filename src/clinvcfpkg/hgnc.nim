import strformat
import strutils
import tables
import logging

from ./utils import logger

type 
  Entrez* = TableRef[string, int]
  Alias* = TableRef[string, string]
  SharedAlias* = TableRef[string, int]

type
  HgncIndex* = ref object
    entrez*: Entrez
    alias*: Alias
    sharedAlias*: SharedAlias

proc newEntrez*(): Entrez =
  result = newTable[string, int]()

proc newAlias*(): Alias =
  result = newTable[string, string]()

proc newSharedAlias*(): SharedAlias =
  result =  newTable[string, int]()

proc newHgncIndex*(): HgncIndex =
  ## Instantiate a new HgncIndex
  ## Each gene that have an Entrez ID is stored in 'entrez' attribute.
  ## Each gene alias is stored in 'alias' attribute as key
  result = HgncIndex(entrez: newEntrez(), alias: newAlias(), sharedAlias: newSharedAlias())

const
  ColApprovedSymbol = "Approved symbol"
  ColAliasSymbols   = "Alias symbols"
  ColEntrezId       = "NCBI Gene ID(supplied by NCBI)"

proc initHgncDbfromFile*(file: string): HgncIndex =
  ## Create an HgncIndex from HGNC table.
  ## The file must contain the required columns but their order does not matter.
  result = newHgncIndex()
  let f = open(file)
  defer: f.close()
  var
    line: string
    cols = initTable[string, int]()  # column name -> index
  # Parse header to build column index map
  if not f.read_line(line):
    raise newException(IOError, "HGNC table is empty")
  for i, name in pairs(line.split("\t")):
    cols[name] = i
  for required in [ColApprovedSymbol, ColAliasSymbols, ColEntrezId]:
    if not cols.hasKey(required):
      raise newException(IOError, fmt"HGNC table is missing required column: {required}")
  let
    iSymbol = cols[ColApprovedSymbol]
    iAlias  = cols[ColAliasSymbols]
    iEntrez = cols[ColEntrezId]
  while f.read_line(line):
    var sl = line.split("\t")
    if sl[iEntrez] == "":
      continue
    result.entrez[sl[iSymbol]] = parseInt(sl[iEntrez])
    if sl[iAlias] != "":
      for g in sl[iAlias].split(", "):
        if result.alias.hasKey(g):
          if not result.sharedAlias.hasKey(g):
            result.sharedAlias[g] = 1
          else:
            inc(result.sharedAlias[g])
        else:
          result.alias[g] = sl[iSymbol]
  # Some aliases are shared between genes — remove them as ambiguous
  for g, c in result.sharedAlias:
    logger.log(lvlInfo, fmt"[initHgncDbfromFile] remove ambiguous shared alias {g}")
    if result.alias.hasKey(g):
      result.alias.del(g)

            


