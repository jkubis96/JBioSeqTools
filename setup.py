from setuptools import setup, find_packages

VERSION = '1.0.1' 
DESCRIPTION = 'JBioSeqTools'
LONG_DESCRIPTION = 'JBioSeqTools is the python library for biological sequence optimization (GC % content & codon frequency) for better expression of different species cells in vivo. It also allows building AAV vectors with the possibility of choosing sequences between ITRs such as transcript, promotor, enhancer, and molecular fluorescence tag. Finally, the user obtains ready for order construct with a whole sequence and visualization. Package description  on https://github.com/jkubis96/JBioSeqTools'

# Setting up
setup(
        name="JBioSeqTools", 
        version=VERSION,
        author="Jakub Kubis",
        author_email="jbiosystem@gmail.com",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=['pandas', 'tqdm', 'matplotlib', 'numpy'],       
        keywords=['sequence', 'optimization', 'vectors', 'AAV', 'GC', 'restriction enzimes'],
        license = 'MIT',
        classifiers = [
            "Development Status :: 3 - Alpha",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: POSIX :: Linux",
        ],
        python_requires='>=3.6',
)

