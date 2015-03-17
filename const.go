package treat

type OrientationType  int8
type BaseCountType    uint64

const FORWARD OrientationType =   1
const REVERSE OrientationType =  -1

const BUCKET_ALIGNMENTS = "alignments"
const BUCKET_TEMPLATES  = "templates"
const BUCKET_SAMPLES    = "samples"
const BUCKET_ESIDX    = "esidx"
