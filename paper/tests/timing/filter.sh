#!/bin/bash

egrep "^time:| of all | insns per cycle " Testing/Temporary/LastTest.log \
  | awk -f filter1.awk | tr -s " " | tr -s "%" | tr "%" "0" \
  | awk -f filter2.awk \
  | awk -f filter3.awk 

