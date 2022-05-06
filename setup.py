from setuptools import setup

setup(
    name='autorate', 
    version='0.1.0', 
    author = 'Eshan King', 
    author_email = '',
    packages=['autorate', "autorate.test", "autorate.test.data"], 
    install_requires = [
      "pandas",
      "pytest",
      "scipy",
      "matplotlib",
      "numpy"
    ],
    include_package_data=True,
    package_data={'': ['data/*.csv']}
)