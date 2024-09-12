using Pkg

function prepare_project()

    Pkg.activate(".") 

    Pkg.instantiate()
    
    packages = ["Symbolics","Glob", "Plots", "DataFrames", "CSV", "CairoMakie", "Dates", "DoubleExponentialFormulas", "FileIO", "Interpolations", "LaTeXStrings", "LinearAlgebra", "PolynomialRoots", "PyCall", "QuadGK","Colors", "StatsBase", "Symbolics"]
    project_deps = keys(Pkg.project().dependencies)
    for pkg in packages
        if !(pkg in project_deps)
            println("Adding $pkg package to the project")
            Pkg.add(pkg)
        else
            println("$pkg is already here in the project :)")
        end
    end

    println("Updating all packages")
    Pkg.update()

    try
        run(`ffmpeg -version`)
        println("ffmpeg is installed.")
    catch
        println("ffmpeg is not installed, installing now...")
        run(`brew update`)
        run(`brew upgrade`)
        run('brew install ffmpeg')
      
    end

    py_packages = ["matplotlib", "Pillow"]
    py_imports = ["matplotlib", "PIL"]

    for (pkg, py_import) in zip(py_packages, py_imports)
        try
            pyimport(py_import)
            println("$pkg is already installed.")
        catch
            println("Installing Python package: $pkg")
            try
                run(`python3 -m pip install $pkg`)
                println("$pkg installed successfully.")
            catch e
                println("Error installing $pkg: $e")
            end
        end
    end
end

# Run the project preparation steps
prepare_project()

println("Project is ready for use.")