using Pkg

ENV["PYTHON"] = ""

Pkg.build("PyCall")

using PyCall

# install trackpy and pandas
packages = ["pandas", "tables", "trackpy"]
run(PyCall.python_cmd(`-m pip install $packages`))