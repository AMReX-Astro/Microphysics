from setuptools import setup, find_packages

setup(name='StarKiller',
      version='0.1',
      description='Python interfaces to StarKiller Microphysics',
      url='https://github.com/starkiller-astro/Microphysics',
      author='The StarKillers',
      author_email='michael.zingale@stonybrook.edu',
      license='BSD',
      packages=find_packages(),
      package_data={"StarKiller": ["burner/*", "eos/*", "network/*", "interfaces/*", "integration/*", "models/*", "examples/*"]},
      install_requires=['numpy', 'matplotlib'],
      zip_safe=False)
