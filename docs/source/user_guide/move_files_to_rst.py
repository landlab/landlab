#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 18:31:09 2019

@author: barnhark
"""

import os
import glob
import subprocess

files = glob.glob("*.md")
for file in files:
    args = ["pandoc",
           "--from=markdown",
           "--to=rst",
           "--output="+file.replace("md", "rst"),
           file]
    subprocess.call(args)
    os.remove(file)
    
files = glob.glob("*.rest")
for file in files:
    os.rename(file, file.replace("rest", "rst"))

