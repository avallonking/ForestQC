from setuptools import setup

config = {
    'description': 'ForestQC: Quality control on genetic variants from next-generation sequencing data using random forest',
    'author': 'Jiajin Li (Albert Lee)',
    'url': 'None',
    'download_url': 'https://github.com/avallonking/ForestQC',
    'author_email': 'lijj36@ucla.edu',
    'version': '1.1.5.7',
    'install_requires': ['pandas','numpy','scikit-learn'],
    'packages': ['ForestQC'],
    'scripts': [],
    'setup_requires': [],
    'entry_points': {'console_scripts': ['ForestQC = ForestQC.__main__:main']},
    'name': 'ForestQC'
}

setup(**config)
