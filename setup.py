import os
from setuptools import setup

# list of all scripts to be included with package
scripts=[]
for s in ['scripts','GIA','OBP','reanalysis','seismic','SMB','TWS']:
    scripts.extend([os.path.join(s,f) for f in os.listdir(s) if f.endswith('.py')])

setup(
    name='model-harmonics',
    scripts=scripts,
)
