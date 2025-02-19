# -*- coding = utf-8 -*-
from setuptools import setup, find_packages # type: ignore

setup(
    name="CircCode3",
    version="1.0",
    author="zzh",
    description="CircCode3 is a software that looks for evidence of targeted circRNA translation from Ribo-seq and mass spectrometry data.",
    url="git url",  # 建议替换为具体的git仓库URL
    license="GPLv3",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPLv3",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(where='src'),
    package_dir={"": "src"},  # 指定 src 目录作为包的根目录
    include_package_data=True,  # 包含 MANIFEST.in 指定的文件
    install_requires=[],  # 如果有依赖项，请在这里列出
    entry_points={
        "console_scripts": [
            "CircCode3=CircCode3.main:main",  # src 目录下的CircCode3.py作为入口点
        ],
    },
)