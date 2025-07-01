from setuptools import (
    setup,
    find_packages
)

from CATS import __version__

setup(
    name="CATS",
    version=__version__,
    packages=find_packages(),
    package_data={"CATS": ["../db/human/*", "../db/mouse/*"]},
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "CATS = CATS.main:main",
            "CATS-converter = gui.converter:main",
        ],
    },
)
