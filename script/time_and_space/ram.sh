#!/usr/bin/bash

while true; do
    python3 -c "import psutil;print(psutil.virtual_memory().percent)"
    sleep 0.5
done

# pip install psutil
# nohup bash ram.sh >> space.csv &