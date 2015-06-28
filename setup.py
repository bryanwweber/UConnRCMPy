from setuptools import setup, find_packages

setup(
    name='UConnRCMPy',
    version='0.1.0',
    description='A package to process RCM data',
    url='https://github.com/bryanwweber/UConnRCMPy',
    author='Bryan W. Weber',
    author_email='weber@engr.uconn.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'process-ignition-loop=uconnrcmpy.ign_loop:main',
        ],
    },
)
