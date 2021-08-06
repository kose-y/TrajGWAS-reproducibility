for effectsize in [0.01, 0.05, 0.1, 0.15, 0.2, 0.25]
    for m in [6000, 100_000] 
        for maf in [0.001, 0.05, 0.3]
            for ni in [2:8, 10:30]
                println("submit job for m=$m, effectsize=$effectsize, maf=$maf, ni=$ni")
                jcode = "effectsize, m, ni, maf = $effectsize, $m, $ni, $maf; include(\"vgwas_power_template_new.jl\")"
                # prepare sh file for qsub
                open("tmp.sh", "w") do io
                    println(io, "#!/bin/bash")
                    println(io, "#\$ -cwd")
                    println(io, "# error = Merged with joblog")
                    println(io, "#\$ -o joblog.\$JOB_ID")
                    println(io, "#\$ -j y")
                    if m == 100_000
                        println(io, "#\$ -l h_rt=6:30:00,h_data=8G,arch=intel*") # request runtime and memory
                    else
                        println(io, "#\$ -l h_rt=4:00:00,h_data=4G,arch=intel*") # request runtime and memory
                    end
                    println(io, "#\$ -pe shared 2") # request runtime and memory
                    println(io, "# Email address to notify")
                    println(io, "#\$ -M \$USER@mail")
                    println(io, "# Notify when")
                    println(io, "#\$ -m a")
                    println(io)
                    println(io, "# load the job environment:")
                    println(io, ". /u/local/Modules/default/init/modules.sh")
                    println(io, "module load julia/1.5.1")
                    println(io)
                    println(io, "# run julia code")
                    println(io, "julia -e '$jcode' > output.\$JOB_ID 2>&1")
                end
                # submit job
                run(`qsub tmp.sh`)
            end
        end
    end
end