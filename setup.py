import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="spktype01",
    version="1.0.0",
    author="whiskie14142",
    author_email="whiskie14142@gmail.com",
    description="A supporting module for jplephem to handle data type 1",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/whiskie14142/spktype01",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)