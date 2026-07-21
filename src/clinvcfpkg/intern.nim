import tables

type
  StringInterner* = object
    ## Deduplicates repetitive strings by mapping each distinct value to a
    ## stable int32 id. Storing the id instead of the string lets millions of
    ## records share a single copy of each value (ORC-safe, unlike string
    ## buffer sharing which arc/orc does not provide).
    toId: Table[string, int32]
    toStr: seq[string]

proc initStringInterner*(): StringInterner =
  result.toId = initTable[string, int32]()
  result.toStr = @[]

proc intern*(si: var StringInterner, s: string): int32 =
  ## Return the id for `s`, assigning a new one the first time it is seen.
  si.toId.withValue(s, id):
    return id[]
  do:
    result = int32(si.toStr.len)
    si.toStr.add(s)
    si.toId[s] = result

proc `[]`*(si: StringInterner, id: int32): string {.inline.} =
  ## Resolve an interned id back to its original string.
  si.toStr[id]