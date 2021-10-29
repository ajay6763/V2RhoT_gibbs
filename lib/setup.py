from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Absolute shear-wave velocity to temperature'
LONG_DESCRIPTION = 'WILL COME!!'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="V2RhoT_gibbs", 
        version=VERSION,
        author="Ajay Kumar",
        author_email="<ajay673@gmail.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=["numpy", "scipy"], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'first package'],
        classifiers= [
            "Development Status :: initial",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 3",
        ]
)