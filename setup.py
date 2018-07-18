from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='esutils',
      version='0.1',
      description='Electronic structure utilities',
      long_description=readme(),
      url='http://github.com/hpmartins/esutils',
      author='H. P. Martins',
      author_email='hpmartins@gmail.com',
      license='MIT',
      packages=['esutils'],
      install_requires=[
          'numpy',
          'scipy',
          'pandas',
      ],
      include_package_data=True,
      zip_safe=False)