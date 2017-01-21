from setuptools import setup, find_packages

setup(
    name='rosemary',
    version='0.0.1',
    url='https://github.com/siegelzero/rosemary',
    license='MIT',
    author='Kenneth Brown',
    author_email='siegel.zero@gmail.com',
    description='Algorithms for combinatorics and number theory',
    long_description="""rosemary is a suite of algorithms for computations in number theory, discrete mathematics, and
                        combinatorial optimization""",
    include_package_data=True,
    # packages=['rosemary', 'rosemary.number_theory'],
    packages=find_packages(),
)
