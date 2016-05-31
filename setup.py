from setuptools import setup, find_packages

setup(
    name='UConnRCMPy',
    version='2.1.0',
    description='A package to process RCM data',
    url='https://github.com/bryanwweber/UConnRCMPy',
    author='Bryan W. Weber',
    author_email='bryan.weber@uconn.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'processrcmfolder=uconnrcmpy.dataprocessing:process_folder',
        ],
    },
)
