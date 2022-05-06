from setuptools import setup

setup(
    name='evodm', 
    version='1.0.0', 
    author = 'Davis Weaver', 
    author_email = 'dtw43@case.edu',
    packages=['autorate', "autorate.test"], 
    install_requires = [
      "pandas",
      "pytest"
    ],
)