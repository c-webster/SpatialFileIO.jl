extension(url::String) = match(r"\.[A-Za-z0-9]+$", url).match

function readlas(infile::String)

    if extension(infile) == ".laz"
        header, dsmdat = LazIO.load(infile)
    elseif extension(infile) == ".las"
        header, dsmdat = FileIO.load(infile)
    else
        error("Unknown DSM file extension")
    end

        dsm_x = fill(NaN,size(dsmdat))
        dsm_y = fill(NaN,size(dsmdat))
        dsm_z = fill(NaN,size(dsmdat))
        dsm_c = fill(NaN,size(dsmdat))

        for dx in eachindex(dsmdat)
            dsm_c[dx] = trunc(Int,float(classification(dsmdat[dx])))
            if dsm_c[dx] == 2
                continue
            else
                dsm_x[dx] = xcoord(dsmdat[dx],header)
                dsm_y[dx] = ycoord(dsmdat[dx],header)
                dsm_z[dx] = zcoord(dsmdat[dx],header)
            end
        end

        rows = findall(isnan,dsm_x)
        deleteat!(dsm_x,rows)
        deleteat!(dsm_y,rows)
        deleteat!(dsm_z,rows)

        return dsm_x, dsm_y, dsm_z

end


function importdtm(dtmf::String,tilt::Bool,limits::Any=nothing)

    if tilt

        file = matopen(dtmf); dtmdat = read(file,"dtm"); close(file)
        dtm_x = dtmdat["x"]
        dtm_y = dtmdat["y"]
        dtm_z = dtmdat["z"]
        dtm_s = dtmdat["s"]
        dtm_a = dtmdat["a"]
        dtm_cellsize = dtmdat["cellsize"]

        return dtm_x, dtm_y, dtm_z, dtm_s, dtm_a, dtm_cellsize

    else

        if extension(dtmf) == ".mat"

            file = matopen(dtmf); dtmdat = read(file,"dtm"); close(file)
            dtm_x = vec(dtmdat["x"])
            dtm_y = vec(dtmdat["y"])
            dtm_z = vec(dtmdat["z"])
            dtm_cellsize = dtmdat["cellsize"]
            rows = findall(isnan,dtm_z)
            deleteat!(dtm_x,rows)
            deleteat!(dtm_y,rows)
            deleteat!(dtm_z,rows)

        elseif extension(dtmf) == ".asc" || extension(dtmf) == ".txt" ||  extension(dtmf) == ".tif"

                dtm_x, dtm_y, dtm_z, dtm_cellsize = read_griddata(dtmf,true,true)

        end

    return dtm_x, dtm_y, dtm_z, dtm_cellsize

    end

end

"""
Read gridded spatial data in .tif or .asc format.
Will also take a file with .txt extention if it is in the same format as an .asc file

Returns x,y,z and cellsize data for grid either as 1D or 2D arrays.

`read_griddata(fname,delete_rows,vectorize)`

# Arguments
- fname::String : full filepath and name of file to be read
- vectorize::Bool : true=return x,y,z data as three single column vectors; false=return x,y,z data as 2D matrices
    default=true
- delete_rows::Bool : to delete rows with NaN values - requires data to be vectorized.
    default=true
    if false, function does not return cellsize.
"""
function read_griddata(fname::String,vectorize=true::Bool,delete_rows=true::Bool)

    if extension(fname) == ".asc" || extension(fname) == ".txt"
        f = open(fname)
            ncols     = parse(Int64,split(readline(f))[2])
            nrows     = parse(Int64,split(readline(f))[2])
            xllcorner = parse(Float64,split(readline(f))[2])
            yllcorner = parse(Float64,split(readline(f))[2])
            cellsize  = parse(Float64,split(readline(f))[2])
            nodatval  = parse(Float64,split(readline(f))[2])
        close(f)

        dat = readdlm(fname,skipstart=6)

        tgrid = Matlab.meshgrid(collect(xllcorner:cellsize:(xllcorner+cellsize*(ncols-1))) .+ cellsize/2,
                                collect(yllcorner:cellsize:(yllcorner+cellsize*(nrows-1))) .+ cellsize/2)

    elseif extension(fname) == ".tif"

        dataset = ArchGDAL.read(fname)

        gt = ArchGDAL.getgeotransform(dataset)
        xulcorner = gt[1]
        cellsize  = gt[2]
        yulcorner = gt[4]
        nodatval  = ArchGDAL.getnodatavalue(ArchGDAL.getband(dataset,1))
        ncols     = ArchGDAL.width(dataset)
        nrows     = ArchGDAL.height(dataset)

        dat = Float64.(transpose(ArchGDAL.read(ArchGDAL.getband(dataset,1))))

        tgrid = Matlab.meshgrid(collect(xulcorner:cellsize:xulcorner+(cellsize*ncols-cellsize)).+cellsize/2,
                                collect(yulcorner-(cellsize*nrows-cellsize):cellsize:yulcorner).+cellsize/2)

    end

    replace!(dat, nodatval=>NaN)

    if vectorize
        dat_x = vec(tgrid[1]);
        dat_y = vec(tgrid[2]);
        dat_z = vec(reverse(dat,dims=1))

        if delete_rows
            rows = findall(isnan,dat_z)
            deleteat!(dat_x,rows)
            deleteat!(dat_y,rows)
            deleteat!(dat_z,rows)
        end

        return dat_x, dat_y, dat_z, cellsize

    else
        return tgrid[1], tgrid[2], dat, cellsize
    end

end

"""
Reads header/meta information for gridded spatial data

# Usage
`ncols, nrows, xllcorner, yllcorner, cellsize, nodatval = read_griddata_header(fname::String)`

If fname is .asc or .txt format, function assumes coordinates are of lower left corner of extent
"""
function read_griddata_header(fname::String)

    if extension(fname) == ".asc" || extension(fname) == ".txt"
        f = open(fname)
            ncols     = parse(Int64,split(readline(f))[2])
            nrows     = parse(Int64,split(readline(f))[2])
            xllcorner = parse(Float64,split(readline(f))[2])
            yllcorner = parse(Float64,split(readline(f))[2])
            cellsize  = parse(Float64,split(readline(f))[2])
            nodatval  = parse(Float64,split(readline(f))[2])
        close(f)

    elseif extension(fname) == ".tif"

        dataset = ArchGDAL.read(fname)

        gt = ArchGDAL.getgeotransform(dataset)
        xllcorner = gt[1] # also xulcorner
        cellsize  = gt[2]
        nodatval  = ArchGDAL.getnodatavalue(ArchGDAL.getband(dataset,1))
        ncols     = ArchGDAL.width(dataset)
        nrows     = ArchGDAL.height(dataset)
        yllcorner = gt[4]-(cellsize*nrows) # gt[4] = yulcorner

    end

    return ncols, nrows, xllcorner, yllcorner, cellsize, nodatval

end


"""
Imports data window from geotiff

# Usage
`dat_x, dat_y, dat_z = read_griddata_window(fname::String,limits::Array{Float64,1},
                                vectorize=true::Bool,delete_rows=true::Bool))`

if vectorize=true, values are returned as 1-D arrays
if delete_rows=true, all NaN values are deleted

Note, if vectorize=false, delete_rows cannot be reached.

limits should be matrix of [xmin xmax ymin ymax]

"""
function read_griddata_window(fname::String,limits,
                                vectorize=true::Bool,delete_rows=true::Bool)

    dataset = ArchGDAL.read(fname)

    gt = ArchGDAL.getgeotransform(dataset)
    xllcorner = gt[1] # also xulcorner
    cellsize  = gt[2]
    nodatval  = ArchGDAL.getnodatavalue(ArchGDAL.getband(dataset,1))
    ncols     = ArchGDAL.width(dataset)
    nrows     = ArchGDAL.height(dataset)
    yllcorner = gt[4]-(cellsize*nrows) # gt[4] = yulcorner


    if ((xllcorner .< limits[1] .< xurcorner) || (xllcorner .< limits[2] .< xurcorner)) &&
            ((yllcorner .< limits[3] .< yurcorner) || (yllcorner .< limits[4] .< yurcorner))

        # check bounds aren't east or south of data
        if limits[1] .< xllcorner
            limits[1] = xllcorner
        end

        if limits[3] .< yllcorner
            limits[3] = yllcorner
        end

        # cell boundaries
        dimsx = xllcorner:cellsize:limits[1]
        dimsy = gt[4]:-cellsize:limits[4]

        # size of window in cells
        xoffset = size(dimsx)[1]-1
        yoffset = size(dimsy)[1]-1

        # correct for negative offset (although should be negated by the above check bounds step)
        if xoffset < 0; xoffset = 0; end
        if yoffset < 0; yoffset = 0; end

        # get size of window in cells
        xmin  = dimsx[end]
        xmax  = (xmin:cellsize:limits[2])[end]
        if xmax - limits[2] !== 0.0; xmax += cellsize; end

        ymin = (yllcorner:cellsize:limits[3])[end]
        ymax = (ymin:cellsize:limits[4])[end]
        if ymax - limits[4] !== 0.0; ymax += cellsize; end

        xsize = (xmax - xmin)/cellsize
        ysize = (ymax - ymin)/cellsize

        # check bounds aren't west or north of data
        if xsize > ncols; xsize = ncols; end
        if ysize > nrows; ysize = nrows; end

        indat = Float64.(transpose(ArchGDAL.read(dataset, 1, Int(xoffset), Int(yoffset), Int(xsize), Int(ysize))))

        tgrid = Matlab.meshgrid(collect(xmin:cellsize:xmax-cellsize) .+ cellsize/2,
                                collect(ymin:cellsize:ymax-cellsize) .+ cellsize/2)

        replace!(indat, nodatval=>NaN)

        if vectorize
            dat_x = vec(tgrid[1]);
            dat_y = vec(tgrid[2]);
            dat_z = vec(reverse(indat,dims=1))

            if delete_rows
                rows = findall(isnan,dat_z)
                deleteat!(dat_x,rows)
                deleteat!(dat_y,rows)
                deleteat!(dat_z,rows)
            end

            return dat_x, dat_y, dat_z, cellsize

        else
            return tgrid[1], tgrid[2], dat, cellsize
        end

    else
        error("Requested window out of bounds of dataset")
    end

end
