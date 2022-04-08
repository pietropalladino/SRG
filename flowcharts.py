from pyflowchart import Flowchart

with open('SRG.py') as f:

   code = f.read()
   fc = Flowchart.from_code(code)
   print(fc.flowchart())
with open('myfunctions.py') as g:
   code=g.read()
   gc=Flowchart.from_code(code)
   print(gc.flowchart())
