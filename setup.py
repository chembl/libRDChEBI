from setuptools import setup

if __name__ == "__main__":
    setup(
        name="libRDChEBI",
        version="0.2.4",
        author="Eloy Félix",
        author_email="chebi-help@ebi.ac.uk",
        description="RDKit library to deal with ChEBI's chemistry",
        url="https://github.com/chembl/libRDChEBI",
        license="MIT",
        packages=["libRDChEBI"],
        long_description=open("README.md", encoding="utf-8").read(),
        long_description_content_type="text/markdown",
        install_requires=["chembl_structure_pipeline~=1.2.0", "rdkit>=2023.09.4"],
        tests_require=["pytest"],
        classifiers=[
            "Intended Audience :: Developers",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
        ],
        zip_safe=True,
    )
