from setuptools import setup, find_packages
 
classifiers = [
  'Development Status :: 5 - Production/Stable',
  'Intended Audience :: Education',
  'Operating System :: Microsoft :: Windows :: Windows 10',
  'License :: OSI Approved :: MIT License',
  'Programming Language :: Python :: 3'
]
 
setup(
  name='mf3cc',
  version='2.0.0',
  description='Model-free 3-component compact polarimetric decomposition',
  py_modules = ['mf3cc'],
  long_description=open('README.txt').read() + '\n\n' + open('CHANGELOG.txt').read(),
  url='',  
  author='Subhadip Dey',
  author_email='sdey2307@gmail.com',
  license='MIT', 
  classifiers=classifiers,
  keywords='SAR, Decomposition, Compact pol, Model-free, MF3CC', 
  packages=find_packages(),
  install_requires=['GDAL', 'numpy'] 
)