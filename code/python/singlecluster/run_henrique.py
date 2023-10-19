import py_compile
py_compile.compile('main.py')
import main

name = 'LaMnO3' 
C = main.run_cluster(name)

name = 'LaFeO3' 
C = main.run_cluster(name)
