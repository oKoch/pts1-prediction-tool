from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='pts1_prediction_tool',
    version='1.0.0',
    description='Tool for predicting a PTS1 (peroxisomal targeting sequence 1) in a aminoacidsequence.',
    license='BSD-3',
    author='Oliver Koch',
    author_email='ollis.lab@gmail.com',
    keywords=['pts1', 'pts1-prediction', 'pts1_prediction_tool'],
    url='https://github.com/oKoch/pts1-prediction-tool',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['pts1_prediction_tool'],
    # package_dir={"": "pts1_prediction_tool"},
    python_requires=">=3.6",
    package_data={
            # If any package contains *.fasta files, include them:
        "": ["*.fasta"],
        # And include any *.dat files found in the "data" subdirectory
        # of the "mypkg" package, also:
        "pts1_prediction_tool": ["datasets/*.fasta"]
    }
)
