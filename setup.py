# Created on Wed Dec 21 15:26:26 2019

# Author: XiaoTao Wang

"""
Setup script for NeoLoopFinder.

"""
import os, sys, neoloop, glob
import setuptools

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if (sys.version_info.major!=3) or (sys.version_info.minor<6):
    print('PYTHON 3.6+ IS REQUIRED. YOU ARE CURRENTLY USING PYTHON {}'.format(sys.version.split()[0]))
    sys.exit(2)

# Guarantee Unix Format
for src in glob.glob('scripts/*'):
    text = open(src, 'r').read().replace('\r\n', '\n')
    open(src, 'w').write(text)

setuptools.setup(
    name = 'neoloop',
    version = neoloop.__version__,
    author = neoloop.__author__,
    author_email = 'wangxiaotao686@gmail.com',
    url = 'https://github.com/XiaoTaoWang/NeoLoopFinder',
    description = 'Predict neo-loops induced by structural variations',
    keywords = 'Hi-C cooler cancer enhancer hijacking',
    long_description = read('README.rst'),
    long_description_content_type='text/x-rst',
    scripts = glob.glob('scripts/*'),
    packages = setuptools.find_packages(),
    package_data = {
        '': ['data/*.pkl']
    },
    classifiers = [
        'Programming Language :: Python :: 3 :: Only',
        'Operating System :: POSIX',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ]
    )
