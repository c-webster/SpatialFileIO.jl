module SpatialFileIO

using LasIO, LazIO, DelimitedFiles, VectorizedRoutines, ArchGDAL

include("FilesIn.jl")
include("FilesOut.jl")

export
    extension,
    readlas,
    importdtm,
    read_ascii,
    read_ascii_header,
    read_geotiff,
    write_ascii

end
