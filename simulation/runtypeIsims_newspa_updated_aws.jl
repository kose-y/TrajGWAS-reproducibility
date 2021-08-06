conditions = [[6000, 6:10, 0.01], 
	      [6000, 6:10, 0.05],
	      [6000, 6:10, 0.30], 
	      [6000, 10:30, 0.01], 
	      [6000, 10:30, 0.05], 
	      [6000, 10:30, 0.30],
              [100_000, 6:10, 0.001], 
	      [100_000, 6:10, 0.05],
	      [100_000, 6:10, 0.30], 
	      [100_000, 10:30, 0.001], 
	      [100_000, 10:30, 0.05], 
	      [100_000, 10:30, 0.30],
	      ] #m, ni, maf
for index in 1:100
    for condition in conditions
        m, ni, maf = condition
        println("submit job for index=$index, maf=$maf, ni=$ni, m=$m")
        jcode = "index, m, ni, maf = $index, $m, $ni, $maf; include(\"vgwastypeIsimtemplate_newspa_updated.jl\")"
        # prepare sh file for qsub
        open("tmp.sh", "w") do io
            println(io, "#!/bin/bash")
            println(io, "#\$ -cwd")
            println(io, "# error = Merged with joblog")
            println(io, "#\$ -o joblog.\$JOB_ID")
            println(io, "#\$ -j y")
            println(io, "#\$ -l h_rt=4:00:00") # request runtime and memory
            println(io, "# Email address to notify")
            println(io, "#\$ -M y3kkoseyoon@gmail.com")
            println(io, "# Notify when")
            println(io, "#\$ -m ea")
            println(io)
            println(io, "# load the job environment:")
            println(io, ". /shared/julia_setup.sh")
            println(io)
            println(io, "export OMP_NUM_THREADS=1")
            println(io, "# run julia code")
            println(io, "julia -e '$jcode' > output.\$JOB_ID 2>&1")
        end
        # submit job
        try
            run(`qsub tmp.sh`)
        catch e
            println(e)
        end
    end
end
