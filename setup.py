from setuptools import find_packages, setup

setup(
    name='ipaPy2',
    packages=find_packages(include=['ipaPy2']),
    version='0.1.5',
    description='Integrated Probabilistic Annotation (IPA) 2.0 - Python implementation ',
    author='Dr Francesco Del Carratore',
    author_email='francescodc87@gmail.com',
    url='https://github.com/francescodc87/ipaPy2',
    license='MIT',
    install_requires=[
          'pandas',
          'molmass==2021.6.18',
          'scipy==1.8.1',
          'tqdm==4.64.0'],
    setup_requires = ['pytest-runner'],
    tests_require= ['pytest==4.4.1'],
    test_suite='test'
)
