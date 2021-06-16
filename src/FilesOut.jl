function write_ascii(headers,fname::String,data,reshape::Bool)

    data[isnan.(data)] .= Int.(headers[6])

    open(fname;write=true) do f
            write(f,"ncols "*string(headers[1])*"\n")
            write(f,"nrows "*string(headers[2])*"\n")
            write(f,"xllcorner "*string(headers[3])*"\n")
            write(f,"yllcorner "*string(headers[4])*"\n")
            write(f,"cellsize "*string(headers[5])*"\n")
            write(f,"NODATA_value "*string(headers[6])*"\n")
    end
    f = open(fname,"a")
            if reshape
                writedlm(f,reverse(reshape(data,Base.parse(Int64,headers[2]),Base.parse(Int64,headers[1])),dims=1))
            else
                writedlm(f,data)
            end
    close(f)

end
