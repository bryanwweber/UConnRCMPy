from setuptools import setup, find_packages
import os.path as op

with open(op.join(op.dirname(op.realpath(__file__)), 'uconnrcmpy', '_version.py')) as version_file:
    exec(version_file.read())

with open(op.join(op.dirname(op.realpath(__file__)), 'README.md')) as readme_file:
    readme = readme_file.read()

with open(op.join(op.dirname(op.realpath(__file__)), 'CHANGELOG.md')) as changelog_file:
    changelog = changelog_file.read()

setup(
    name='UConnRCMPy',
    version=__version__,
    description='A package to process RCM data',
    long_description=readme + '\n\n' + changelog,
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
        'Programming Language :: Python :: 3.6',
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'processrcmfolder=uconnrcmpy.dataprocessing:process_folder',
        ],
    },
)
