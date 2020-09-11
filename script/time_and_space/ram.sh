#!/usr/bin/bash
# Caution: compiler varies on pc 

while true; do
    python3 -c "import psutil;print(psutil.virtual_memory().percent)"
    sleep 0.5
done

