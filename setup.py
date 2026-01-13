"""

SPRINT

"""

from setuptools import setup, find_packages
import os

# 读取README.md作为长描述
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="sprint",
    version="1.0.0",
    author="Minghan Li",
    author_email="limh25@m.fudan.edu.cn",
    description="A toolkit for accurately identifying RNA editing sites without the need to filter SNPs based on python3",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-username/your-repository-name",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "requests>=2.31.0",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # 与你的LICENSE文件对应
        "Operating System :: OS Independent",
    ],
    # 关键：添加这部分配置，生成sprint命令行入口
    entry_points={
        'console_scripts': [
            'sprint = sprint.__init__:main',  
            'sprint_prepare = sprint.sprint_prepare:main', 
            'sprint_main = sprint.sprint_main:main', 
            'sprint_main_parallel = sprint.sprint_main_parallel:main', 
            'sprint_from_bam = sprint.sprint_from_bam:main', 
            'sprint_from_bam_parallel = sprint.sprint_from_bam_parallel:main'

        ],
    }
)