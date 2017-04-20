from setuptools import setup, find_packages
import os.path as op
import sys

with open(op.join(op.dirname(op.realpath(__file__)), 'uconnrcmpy', '_version.py')) as version_file:
    exec(version_file.read())

with open(op.join(op.dirname(op.realpath(__file__)), 'README.md')) as readme_file:
    readme = readme_file.read()

with open(op.join(op.dirname(op.realpath(__file__)), 'CHANGELOG.md')) as changelog_file:
    changelog = changelog_file.read()

with open(op.join(op.dirname(op.realpath(__file__)), 'CITATION.md')) as citation_file:
    citation = citation_file.read()

desc = readme + '\n\n' + changelog + '\n\n' + citation
try:
    import pypandoc
    long_description = pypandoc.convert_text(desc, 'rst', format='md')
    with open(op.join(op.dirname(op.realpath(__file__)), 'README.rst'), 'w') as rst_readme:
        rst_readme.write(long_description)
except (ImportError, OSError, IOError):
    long_description = desc

install_requires = [
    'cantera>=2.3.0',
    'numpy>=1.8.0',
    'scipy>=0.18.0',
    'pyyaml>-3.12',
    'matplotlib>=1.4.0',
    'pyperclip>=1.5.27'
]

tests_require = [
    'pytest>=3.0.0',
    'pytest-cov>=2.3.1',
]

needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
setup_requires = ['pytest-runner'] if needs_pytest else []

setup(
    name='UConnRCMPy',
    version=__version__,
    description='A package to process RCM data',
    long_description=long_description,
    url='https://github.com/bryanwweber/UConnRCMPy',
    author='Bryan W. Weber',
    author_email='bryan.weber@uconn.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'processrcmfolder=uconnrcmpy.dataprocessing:process_folder',
        ],
    },
    install_requires=install_requires,
    tests_require=tests_require,
    setup_requires=setup_requires,
)
