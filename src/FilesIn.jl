extension(url::String) = match(r"\.[A-Za-z0-9]+$", url).match


"""
Imports data window from las or laz file

# Usage
`dat_x, dat_y, dat_z = readlas(fname::String,limits::Any=nothing))`

# Arguments
- fname::String : full filepath and name of file to be read
- limits::Matrix : Matrix of [xmin xmax ymin ymax] describing boundary of window to load

Notes:
 - currently only loads lidar points classed as 3, 4 or 5 (vegetation classes)

"""
function readlas(fname::String,limits::Any=nothing)

	if extension(fname) == ".laz"
		header, lasdat = LazIO.load(fname)
	elseif extension(fname) == ".las"
		header, lasdat = FileIO.load(fname)
	else
		error("Unknown las file extension")
	end

	if limits == nothing
		limits = [header.x_min,header.x_max,header.y_min,header.y_max]
	end

	dat = DataFrame(lasdat)

	las_c = Int.(dat.raw_classification)

	dx = ((limits[1] .<= (dat.x .* header.x_scale .+ header.x_offset) .<= limits[2]) .&
			(limits[3] .<= (dat.y .* header.y_scale .+ header.y_offset) .<= limits[4])) .&
	 		((las_c .== 3) .| (las_c .== 4) .| (las_c .== 5))

	las_x = dat.x[dx] .* header.x_scale .+ header.x_offset
	las_y = dat.y[dx] .* header.y_scale .+ header.y_offset
	las_z = dat.z[dx] .* header.z_scale .+ header.z_offset

	rows = findall(isnan,las_x)
	deleteat!(las_x,rows)
	deleteat!(las_y,rows)
	deleteat!(las_z,rows)

	return las_x, las_y, las_z

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

        tgrid = (((xllcorner:cellsize:(xllcorner+cellsize*(ncols-1)))' .* ones(nrows)) .+ cellsize/2),
                        ((ones(ncols)' .* (yllcorner:cellsize:(yllcorner+cellsize*(nrows-1)))) .+ cellsize/2)

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

        tgrid = (((xulcorner:cellsize:xulcorner+(cellsize*ncols-cellsize))'  .* ones(nrows)) .+ cellsize/2),
                    ((ones(ncols)' .* (yulcorner-(cellsize*nrows-cellsize):cellsize:yulcorner)) .+ cellsize/2)

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
    if false, function does not return cellsize.
- remove_zeros::Bool : removes all null values in data (useful for normalised data)
    default=false

Notes:
 - if vectorize=false, delete_rows cannot be reached (data stays in 2D).
 - care should be taken if using datasets with cellsize precision greater than 2 d.p.

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

        xsize = Int((xmax - xmin)/cellsize)
        ysize = Int((ymax - ymin)/cellsize)

        # check bounds aren't west or north of data
        if xsize > ncols; xsize = ncols; end
        if ysize > nrows; ysize = nrows; end

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

    else
        error("SpatialFileIO bounds error: Requested window outside bounds of dataset")
    end

end
