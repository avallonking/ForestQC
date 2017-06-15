try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Random forest classification on good/bad variants from vcf files',
    'author': 'Albert Lee',
    'url': 'None',
    'download_url': 'None',
    'author_email': 'lijj36@ucla.edu',
    'version': 'alpha',
    'install_requires': ['nose','pandas','numpy','sklearn'],
    'packages': ['machine_learning_project'],
    'scripts': [],
    'name': 'variants_classifier'
}

setup(**config)