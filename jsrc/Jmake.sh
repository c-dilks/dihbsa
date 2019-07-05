#!/bin/bash

if [ -z "$CLASSPATH" ]; then
  echo "error: CLASSPATH not set; you should source env.sh"
  exit
fi


echo "[+] compiling jsrc..."
javac simpleAnalyze.java

echo "[+] done"
