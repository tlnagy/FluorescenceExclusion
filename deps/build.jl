using PyCall

# install trackpy and pandas
packages = ["pandas", "tables", "trackpy"]
run(PyCall.python_cmd(`-m pip install --user $packages`))