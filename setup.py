from setuptools import setup

setup(
    name='autorate', 
    version='0.1.0', 
    author = 'Eshan King', 
    author_email = '',
    packages=['autorate', "autorate.test"], 
    install_requires = [
      "pandas",
      "pytest"
    ],
)