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


function importdtm(dtmf::String,tilt::Bool)
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
        elseif extension(dtmf) == ".asc" || extension(dtmf) == ".txt"
                dtm_x, dtm_y, dtm_z, dtm_cellsize = read_ascii(dtmf,true)
        end
    return dtm_x, dtm_y, dtm_z, dtm_cellsize
    end
end



function read_ascii(fname::String,delete_rows=true::Bool)

        f = open(fname)
        ncols     = parse(Int64,split(readline(f))[2])
        nrows     = parse(Int64,split(readline(f))[2])
        xllcorner = parse(Float64,split(readline(f))[2])
        yllcorner = parse(Float64,split(readline(f))[2])
        cellsize  = parse(Float64,split(readline(f))[2])
        nodatval  = parse(Float64,split(readline(f))[2])
        close(f)

        dat = readdlm(fname,skipstart=6)
        replace!(dat, -9999=>NaN)

        # GC.gc()
        xdat = collect(xllcorner:cellsize:(xllcorner+cellsize*(ncols-1))) .+ cellsize/2;
        ydat = collect(yllcorner:cellsize:(yllcorner+cellsize*(nrows-1))) .+ cellsize/2;
        tgrid = Matlab.meshgrid(xdat,ydat)

        dat_x = vec(tgrid[1]);
        dat_y = vec(tgrid[2]);
        dat_z = vec(reverse(dat,dims=1));

        if delete_rows
            rows = findall(isnan,dat_z)
            deleteat!(dat_x,rows)
            deleteat!(dat_y,rows)
            deleteat!(dat_z,rows)
        end

    return dat_x, dat_y, dat_z, cellsize

end
