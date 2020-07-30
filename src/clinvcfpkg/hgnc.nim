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

proc initHgncDbfromFile*(file: string): HgncIndex =
  ## Create an HgncIndex from HGNC table
  result = newHgncIndex()
  let f = open(file)
  defer: f.close()
  var line : string
  var isHeader = false
  while f.read_line(line):
    if line == "Approved symbol\tAlias symbols\tNCBI Gene ID(supplied by NCBI)":
      isHeader = true
      continue
    if not isHeader:
      raise newException(IOError, "wrong HGNC table header")
    var sl = line.split("\t")
    if sl[2] == "":
      # next if no entrezID is defined
      continue
    result.entrez[sl[0]] = parseInt(sl[2])
    # handle alias
    if sl[1] != "":
      var sAlias = sl[1].split(", ")
      for g in sAlias:
        if result.alias.hasKey(g):
          # store these alias
          if not result.sharedAlias.hasKey(g):
            result.sharedAlias[g] = 1
          else:
            inc(result.sharedAlias[g])
        else:
          result.alias[g] = sl[0]
  # some alias are shared between genes. These alias are ambiguous so we remove them
  for g, c in result.sharedAlias:
    logger.log(lvlInfo, fmt"[initHgncDbfromFile] remove ambiguous shared alias {g}")
    if result.alias.hasKey(g):
      result.alias.del(g)

            


