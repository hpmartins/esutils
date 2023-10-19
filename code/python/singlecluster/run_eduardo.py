import py_compile
py_compile.compile('main.py')
import main

name = 'Rh2O3'


C = main.run_cluster(name)
