from setuptools import setup, find_packages

setup(
    name='phytochempy',
    version='0.1',
    packages=find_packages(),
    package_data={},
    install_requires=[
        'pandas>=2.1.4',
        'numpy>=1.26'
    ],
    extras_require={
    },
    url='https://github.com/alrichardbollans/PhytoChemicalDiversity',
    license='GNU v.3',
    author='Adam Richard-Bollans',
    description='A package for analysing chemical diversity in vascular plants',
    long_description=open('README.md', encoding="utf8").read()
)
