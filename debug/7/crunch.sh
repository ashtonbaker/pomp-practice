#!/bin/bash

clear;

R CMD BATCH script.R ./output/output.Rout;

ls -lh ./output >> ./output/message.txt;

python finished.py
