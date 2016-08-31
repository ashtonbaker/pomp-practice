#!/bin/bash

clear;

R CMD BATCH abridg.R ./output/output.Rout;

ls -lh ./output >> ./output/message.txt;

python finished.py
