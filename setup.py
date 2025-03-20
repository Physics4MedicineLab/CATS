from setuptools import (
    setup,
    find_packages
)

setup(
    name="CATS",
    version="1.0.0",
    packages=find_packages(),
    package_data={"CATS": ["../db/human/*", "../db/mouse/*"]},
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "CATS = CATS.main:main",
        ],
    },
)
