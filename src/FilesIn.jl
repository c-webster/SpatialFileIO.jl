
"""
Imports data from las or laz file with option to limit returned data by area or classification

# Usage
`dat_x, dat_y, dat_z = readlas(fname::String,limits::Any=nothing,keep_ground::Bool=false,keep_all::=false)`

# Arguments
- fname::String : full filepath and name of file to be read
- limits::Matrix : Matrix of [xmin xmax ymin ymax] describing boundary of window to load
- keep_ground = true; loads all points classified 2-5
- keep_all    = true; loads all points, regardless of classification

Notes:
- keep_ground and keep_all cannot both be true. if both are true, keep_ground takes precedent. 
- if keep_ground and keep_all are both false, function loads only canopy points i.e. points classed as 3, 4 or 5 (vegetation classes)

"""
function readlas(fname::String,limits::Any=nothing,keep_ground::Bool=false,keep_all::Bool=false)

    pc = getfield(PointCloud(fname; attributes = (classification)), :data)

    if keep_ground
        pointmask = ((pc.classification .== 2) .| (pc.classification .== 3) .| 
                    (pc.classification .== 4) .| (pc.classification .== 5))
    elseif keep_all
        pointmask = Bool.(ones(size(pc.classification)))
    else
        pointmask = ((pc.classification .== 3) .| (pc.classification .== 4) .| 
                    (pc.classification .== 5))
    end

    if limits == nothing
        coordmask = Bool.(ones(size(pc.classification)))
    else
        coordmask = (pc.x .>= limits[1]) .& (pc.x .<= limits[2]) .&
                (pc.y .>= limits[3]) .& (pc.y .<= limits[4])
    end

    keepat!(pc, pointmask.*coordmask)

	return pc.x, pc.y, pc.z

end


"""
Read gridded spatial data in .tif or .asc format.
Will also take a file with .txt extention if it is in the same format as an .asc file

Returns x,y,z and cellsize data for grid either as 1D or 2D arrays.
x,y values are the centre of the grid cell

`read_griddata(fname,vectorize,delete_rows,remove_zeros)`

# Arguments
- fname::String : full filepath and name of file to be read
- vectorize::Bool : true=return x,y,z data as three single column vectors; false=return x,y,z data as 2D matrices
    default=true
- delete_rows::Bool : to delete rows with NaN values - requires data to be vectorized.
    default=true
    if false, function does not return cellsize.
- remove_zeros::Bool : removes all null values in data (useful for normalised data)
    default=false

"""
function read_griddata(fname::String,vectorize=true::Bool,
                        delete_rows=true::Bool,remove_zeros=false::Bool)

    if splitext(fname)[2] == ".asc" || splitext(fname)[2] == ".txt"
        f = open(fname)
            ncols     = parse(Int64,split(readline(f))[2])
            nrows     = parse(Int64,split(readline(f))[2])
            xllcorner = parse(Float64,split(readline(f))[2])
            yllcorner = parse(Float64,split(readline(f))[2])
            cellsize  = parse(Float64,split(readline(f))[2])
            nodatval  = parse(Float64,split(readline(f))[2])
        close(f)

        dat = readdlm(fname,skipstart=6)

        tgrid = (((xllcorner:cellsize:(xllcorner+cellsize*(ncols-1)))' .* ones(nrows)) .+ cellsize/2),
                        ((ones(ncols)' .* (yllcorner:cellsize:(yllcorner+cellsize*(nrows-1)))) .+ cellsize/2)

    elseif splitext(fname)[2] == ".tif" || splitext(fname)[2] == ".TIFF"

        dataset = ArchGDAL.read(fname)

        gt = ArchGDAL.getgeotransform(dataset)
        xulcorner = gt[1]
        cellsize  = gt[2]
        yulcorner = gt[4]
        nodatval  = ArchGDAL.getnodatavalue(ArchGDAL.getband(dataset,1))
        ncols     = ArchGDAL.width(dataset)
        nrows     = ArchGDAL.height(dataset)

        dat = Float64.(transpose(ArchGDAL.read(ArchGDAL.getband(dataset,1))))

        tgrid = (((xulcorner:cellsize:xulcorner+(cellsize*ncols-cellsize))'  .* ones(nrows)) .+ cellsize/2),
                    ((ones(ncols)' .* (yulcorner-(cellsize*nrows-cellsize):cellsize:yulcorner)) .+ cellsize/2)

    else

        error("File extension not recognised")

    end

    if remove_zeros
		replace!(dat,0=>NaN)
	end

    replace!(dat, nodatval=>NaN)

    if vectorize
        dat_x = vec(tgrid[1]);
        dat_y = vec(tgrid[2]);
        dat_z = vec(reverse(dat,dims=1)) # reversed because the y-axis is inverted

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

Notes:
 - care should be taken if using datasets with cellsize precision greater than 2 d.p.

"""
function read_griddata_header(fname::String)

    if splitext(fname)[2] == ".asc" || splitext(fname)[2] == ".txt"
        f = open(fname)
            ncols     = parse(Int64,split(readline(f))[2])
            nrows     = parse(Int64,split(readline(f))[2])
            xllcorner = parse(Float64,split(readline(f))[2])
            yllcorner = parse(Float64,split(readline(f))[2])
            cellsize  = parse(Float64,split(readline(f))[2])
            nodatval  = parse(Float64,split(readline(f))[2])
        close(f)

    elseif splitext(fname)[2] == ".tif"

        dataset = ArchGDAL.read(fname)

        gt = ArchGDAL.getgeotransform(dataset)
        xllcorner = gt[1] # also xulcorner
        cellsize  = round(gt[2],digits=2)
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
`dat_x, dat_y, dat_z = read_griddata_window(fname::String,limits,vectorize=true::Bool,
                                                delete_rows=true::Bool,remove_zeros=false::Bool)`

# Arguments
- fname::String : full filepath and name of file to be read
- limits::Array : Array of [xmin xmax ymin ymax] describing boundary of window to load
- vectorize::Bool : true=return x,y,z data as three single column vectors; false=return x,y,z data as 2D matrices
    default=true
- delete_rows::Bool : to delete rows with NaN values - requires data to be vectorized.
    default=true
- remove_zeros::Bool : removes all zero values in data and replaces with nan (useful for normalised data)
    default=false

Notes:
 - if vectorize=false, delete_rows cannot be reached (data stays in 2D).
 - care should be taken if using datasets with cellsize precision greater than 2 d.p.
 - if the desired window falls on the edge of the dataset, the function will issue a warning 

"""
function read_griddata_window(fname::String,limits,vectorize=true::Bool,
                                delete_rows=true::Bool,remove_zeros=false::Bool)

    dataset = ArchGDAL.read(fname)

    gt = ArchGDAL.getgeotransform(dataset)
    xllcorner = gt[1] # also xulcorner
    cellsize  = round(gt[2],digits=2)
    nodatval  = ArchGDAL.getnodatavalue(ArchGDAL.getband(dataset,1))
    ncols     = Int(ArchGDAL.width(dataset))
    nrows     = Int(ArchGDAL.height(dataset))
    yllcorner = gt[4]-(cellsize*nrows) # gt[4] = yulcorner
    xurcorner = xllcorner+(cellsize*ncols)
    yurcorner = gt[4]

    # check that at least one corner of the requested limits is within the dataset bounds
    corners = [
        (limits[1], limits[3]),  # bottom-left
        (limits[2], limits[3]),  # bottom-right
        (limits[1], limits[4]),  # top-left
        (limits[2], limits[4])   # top-right
    ]

    within_bounds = any(corners) do (x, y)
        xllcorner <= x <= xurcorner && yllcorner <= y <= yurcorner
    end

    if !within_bounds
        error("SpatialFileIO bounds error: Requested window completely outside bounds of dataset")
    end

    if (!(xllcorner .< limits[1] .< xurcorner) || !(xllcorner .< limits[2] .< xurcorner)) ||
            (!(yllcorner .< limits[3] .< yurcorner) || !(yllcorner .< limits[4] .< yurcorner))

        @warn("Warning: Specified limits extend beyond data bounds. Data will be loaded to available extent.")

        # check individual bounds aren't east or south of data
        (limits[1] .< xllcorner) && (limits[1] = xllcorner)
        (limits[2] .> xurcorner) && (limits[2] = xurcorner)
        (limits[3] .< yllcorner) && (limits[3] = yllcorner)
        (limits[4] .> yurcorner) && (limits[4] = yurcorner)

    end

    # cell boundaries
    dimsx = xllcorner:cellsize:limits[1]
    dimsy = yurcorner:-cellsize:limits[4]

    # offset of window from upper left corner
    xoffset = Int(size(dimsx)[1]-1)
    yoffset = Int(size(dimsy)[1]-1)

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

    # get size of window in cells
    xsize = Int((xmax - xmin)/cellsize)
    ysize = Int((ymax - ymin)/cellsize)

    indat = Float64.(transpose(ArchGDAL.read(dataset, 1, xoffset, yoffset, xsize, ysize)))

    tgrid  = ((xmin:cellsize:xmax-cellsize)'  .* ones(ysize)) .+ cellsize/2,
                (ones(xsize)' .* (ymin:cellsize:ymax-cellsize)) .+ cellsize/2;

    if remove_zeros
        replace!(indat,0=>NaN)
    end

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
        return tgrid[1], tgrid[2], indat, cellsize
    end
    
end
